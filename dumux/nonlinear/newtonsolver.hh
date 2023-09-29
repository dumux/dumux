// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Nonlinear
 * \brief Reference implementation of a Newton solver.
 */
#ifndef DUMUX_NEWTON_SOLVER_HH
#define DUMUX_NEWTON_SOLVER_HH

#include <cmath>
#include <memory>
#include <iostream>
#include <type_traits>
#include <algorithm>
#include <numeric>

#include <dune/common/timer.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpicommunication.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/std/type_traits.hh>
#include <dune/common/indices.hh>
#include <dune/common/hybridutilities.hh>

#include <dune/istl/bvector.hh>
#include <dune/istl/multitypeblockvector.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/typetraits/vector.hh>
#include <dumux/common/typetraits/isvalid.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/common/pdesolver.hh>
#include <dumux/common/variablesbackend.hh>

#include <dumux/io/format.hh>

#include <dumux/linear/matrixconverter.hh>
#include <dumux/assembly/partialreassembler.hh>

#include "newtonconvergencewriter.hh"
#include "primaryvariableswitchadapter.hh"

namespace Dumux::Detail::Newton {

// Helper boolean that states if the assembler exports grid variables
template<class Assembler> using AssemblerGridVariablesType = typename Assembler::GridVariables;
template<class Assembler>
inline constexpr bool assemblerExportsGridVariables
    = Dune::Std::is_detected_v<AssemblerGridVariablesType, Assembler>;

// helper struct to define the variables on which the privarswitch should operate
template<class Assembler, bool exportsGridVars = assemblerExportsGridVariables<Assembler>>
struct PriVarSwitchVariablesType { using Type = typename Assembler::GridVariables; };

// if assembler does not export them, use an empty class. These situations either mean
// that there is no privarswitch, or, it is handled by a derived implementation.
template<class Assembler>
struct PriVarSwitchVariablesType<Assembler, false>
{ using Type = struct EmptyGridVariables {}; };

// Helper alias to deduce the variables types used in the privarswitch adapter
template<class Assembler>
using PriVarSwitchVariables
    = typename PriVarSwitchVariablesType<Assembler, assemblerExportsGridVariables<Assembler>>::Type;

//! helper struct detecting if an assembler supports partial reassembly
struct supportsPartialReassembly
{
    template<class Assembler>
    auto operator()(Assembler&& a)
    -> decltype(a.assembleJacobianAndResidual(std::declval<const typename Assembler::SolutionVector&>(),
                                              std::declval<const PartialReassembler<Assembler>*>()))
    {}
};

// helpers to implement max relative shift
template<class C> using dynamicIndexAccess = decltype(std::declval<C>()[0]);
template<class C> using staticIndexAccess = decltype(std::declval<C>()[Dune::Indices::_0]);
template<class C> static constexpr auto hasDynamicIndexAccess = Dune::Std::is_detected<dynamicIndexAccess, C>{};
template<class C> static constexpr auto hasStaticIndexAccess = Dune::Std::is_detected<staticIndexAccess, C>{};

template<class V, class Scalar, class Reduce, class Transform>
auto hybridInnerProduct(const V& v1, const V& v2, Scalar init, Reduce&& r, Transform&& t)
-> std::enable_if_t<hasDynamicIndexAccess<V>(), Scalar>
{
    return std::inner_product(v1.begin(), v1.end(), v2.begin(), init, std::forward<Reduce>(r), std::forward<Transform>(t));
}

template<class V, class Scalar, class Reduce, class Transform>
auto hybridInnerProduct(const V& v1, const V& v2, Scalar init, Reduce&& r, Transform&& t)
-> std::enable_if_t<hasStaticIndexAccess<V>() && !hasDynamicIndexAccess<V>(), Scalar>
{
    using namespace Dune::Hybrid;
    forEach(std::make_index_sequence<V::N()>{}, [&](auto i){
        init = r(init, hybridInnerProduct(v1[Dune::index_constant<i>{}], v2[Dune::index_constant<i>{}], init, std::forward<Reduce>(r), std::forward<Transform>(t)));
    });
    return init;
}

// Maximum relative shift at a degree of freedom.
// For (primary variables) values below 1.0 we use
// an absolute shift.
template<class Scalar, class V>
auto maxRelativeShift(const V& v1, const V& v2)
-> std::enable_if_t<Dune::IsNumber<V>::value, Scalar>
{
    using std::abs; using std::max;
    return abs(v1 - v2)/max<Scalar>(1.0, abs(v1 + v2)*0.5);
}

// Maximum relative shift for generic vector types.
// Recursively calls maxRelativeShift until Dune::IsNumber is true.
template<class Scalar, class V>
auto maxRelativeShift(const V& v1, const V& v2)
-> std::enable_if_t<!Dune::IsNumber<V>::value, Scalar>
{
    return hybridInnerProduct(v1, v2, Scalar(0.0),
        [](const auto& a, const auto& b){ using std::max; return max(a, b); },
        [](const auto& a, const auto& b){ return maxRelativeShift<Scalar>(a, b); }
    );
}

template<class To, class From>
void assign(To& to, const From& from)
{
    if constexpr (hasStaticIndexAccess<To>() && hasStaticIndexAccess<To>() && !hasDynamicIndexAccess<From>() && !hasDynamicIndexAccess<From>())
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<To::N()>{}, [&](auto i){
            assign(to[Dune::index_constant<i>{}], from[Dune::index_constant<i>{}]);
        });
    }

    else if constexpr (std::is_assignable<To&, From>::value)
        to = from;

    else if constexpr (hasDynamicIndexAccess<To>() && hasDynamicIndexAccess<From>())
        for (decltype(to.size()) i = 0; i < to.size(); ++i)
            assign(to[i], from[i]);

    else if constexpr (hasDynamicIndexAccess<To>() && Dune::IsNumber<From>::value)
    {
        assert(to.size() == 1);
        assign(to[0], from);
    }

    else if constexpr (Dune::IsNumber<To>::value && hasDynamicIndexAccess<From>())
    {
        assert(from.size() == 1);
        assign(to, from[0]);
    }

    else
        DUNE_THROW(Dune::Exception, "Values are not assignable to each other!");
}

} // end namespace Dumux::Detail::Newton

namespace Dumux {

/*!
 * \ingroup Nonlinear
 * \brief An implementation of a Newton solver
 * \tparam Assembler the assembler
 * \tparam LinearSolver the linear solver
 * \tparam Comm the communication object used to communicate with all processes
 * \note If you want to specialize only some methods but are happy with the
 *       defaults of the reference solver, derive your solver from
 *       this class and simply overload the required methods.
 */
template <class Assembler, class LinearSolver,
          class Reassembler = PartialReassembler<Assembler>,
          class Comm = Dune::Communication<Dune::MPIHelper::MPICommunicator> >
class NewtonSolver : public PDESolver<Assembler, LinearSolver>
{
    using ParentType = PDESolver<Assembler, LinearSolver>;

protected:
    using Backend = VariablesBackend<typename ParentType::Variables>;
    using SolutionVector = typename Backend::DofVector;
    using ResidualVector = typename Assembler::ResidualType;
    using LinearAlgebraNativeBackend = VariablesBackend<ResidualVector>;
private:
    using Scalar = typename Assembler::Scalar;
    using JacobianMatrix = typename Assembler::JacobianMatrix;
    using ConvergenceWriter = ConvergenceWriterInterface<SolutionVector, ResidualVector>;
    using TimeLoop = TimeLoopBase<Scalar>;

    // enable models with primary variable switch
    // TODO: Always use ParentType::Variables once we require assemblers to export variables
    static constexpr bool assemblerExportsVariables = Detail::PDESolver::assemblerExportsVariables<Assembler>;
    using PriVarSwitchVariables
        = std::conditional_t<assemblerExportsVariables,
                             typename ParentType::Variables,
                             Detail::Newton::PriVarSwitchVariables<Assembler>>;
    using PrimaryVariableSwitchAdapter = Dumux::PrimaryVariableSwitchAdapter<PriVarSwitchVariables>;

public:
    using typename ParentType::Variables;
    using Communication = Comm;

    /*!
     * \brief The Constructor
     */
    NewtonSolver(std::shared_ptr<Assembler> assembler,
                 std::shared_ptr<LinearSolver> linearSolver,
                 const Communication& comm = Dune::MPIHelper::getCommunication(),
                 const std::string& paramGroup = "")
    : ParentType(assembler, linearSolver)
    , endIterMsgStream_(std::ostringstream::out)
    , comm_(comm)
    , paramGroup_(paramGroup)
    , priVarSwitchAdapter_(std::make_unique<PrimaryVariableSwitchAdapter>(paramGroup))
    {
        initParams_(paramGroup);

        // set the linear system (matrix & residual) in the assembler
        this->assembler().setLinearSystem();

        // set a different default for the linear solver residual reduction
        // within the Newton the linear solver doesn't need to solve too exact
        this->linearSolver().setResidualReduction(getParamFromGroup<Scalar>(paramGroup, "LinearSolver.ResidualReduction", 1e-6));

        // initialize the partial reassembler
        if (enablePartialReassembly_)
            partialReassembler_ = std::make_unique<Reassembler>(this->assembler());
    }

    //! the communicator for parallel runs
    const Communication& comm() const
    { return comm_; }

    /*!
     * \brief Set the maximum acceptable difference of any primary variable
     * between two iterations for declaring convergence.
     *
     * \param tolerance The maximum relative shift between two Newton
     *                  iterations at which the scheme is considered finished
     */
    void setMaxRelativeShift(Scalar tolerance)
    { shiftTolerance_ = tolerance; }

    /*!
     * \brief Set the maximum acceptable absolute residual for declaring convergence.
     *
     * \param tolerance The maximum absolute residual at which
     *                  the scheme is considered finished
     */
    void setMaxAbsoluteResidual(Scalar tolerance)
    { residualTolerance_ = tolerance; }

    /*!
     * \brief Set the maximum acceptable residual norm reduction.
     *
     * \param tolerance The maximum reduction of the residual norm
     *                  at which the scheme is considered finished
     */
    void setResidualReduction(Scalar tolerance)
    { reductionTolerance_ = tolerance; }

    /*!
     * \brief Set the number of iterations at which the Newton method
     *        should aim at.
     *
     * This is used to control the time-step size. The heuristic used
     * is to scale the last time-step size by the deviation of the
     * number of iterations used from the target steps.
     *
     * \param targetSteps Number of iterations which are considered "optimal"
     */
    void setTargetSteps(int targetSteps)
    { targetSteps_ = targetSteps; }

    /*!
     * \brief Set the number of minimum iterations for the Newton
     *        method.
     *
     * \param minSteps Minimum number of iterations
     */
    void setMinSteps(int minSteps)
    { minSteps_ = minSteps; }

    /*!
     * \brief Set the number of iterations after which the Newton
     *        method gives up.
     *
     * \param maxSteps Number of iterations after we give up
     */
    void setMaxSteps(int maxSteps)
    { maxSteps_ = maxSteps; }

    /*!
     * \brief Run the Newton method to solve a non-linear system.
     *        Does time step control when the Newton fails to converge
     * \param vars The variables object representing the current state of the
     *             numerical solution (primary and possibly secondary variables).
     * \param timeLoop The time loop.
     */
    void solve(Variables& vars, TimeLoop& timeLoop) override
    {
        if constexpr (!assemblerExportsVariables)
        {
            if (this->assembler().isStationaryProblem())
                DUNE_THROW(Dune::InvalidStateException, "Using time step control with stationary problem makes no sense!");
        }

        // try solving the non-linear system
        for (std::size_t i = 0; i <= maxTimeStepDivisions_; ++i)
        {
            // linearize & solve
            const bool converged = solve_(vars);

            if (converged)
                return;

            else if (!converged && i < maxTimeStepDivisions_)
            {
                if constexpr (assemblerExportsVariables)
                    DUNE_THROW(Dune::NotImplemented, "Time step reset for new assembly methods");
                else
                {
                    // set solution to previous solution & reset time step
                    Backend::update(vars, this->assembler().prevSol());
                    this->assembler().resetTimeStep(Backend::dofs(vars));

                    if (verbosity_ >= 1)
                    {
                        const auto dt = timeLoop.timeStepSize();
                        std::cout << Fmt::format("Newton solver did not converge with dt = {} seconds. ", dt)
                                  << Fmt::format("Retrying with time step of dt = {} seconds.\n", dt*retryTimeStepReductionFactor_);
                    }

                    // try again with dt = dt * retryTimeStepReductionFactor_
                    timeLoop.setTimeStepSize(timeLoop.timeStepSize() * retryTimeStepReductionFactor_);
                }
            }

            else
            {
                DUNE_THROW(NumericalProblem,
                    Fmt::format("Newton solver didn't converge after {} time-step divisions; dt = {}.\n",
                                maxTimeStepDivisions_, timeLoop.timeStepSize()));
            }
        }
    }

    /*!
     * \brief Run the Newton method to solve a non-linear system.
     *        The solver is responsible for all the strategic decisions.
     * \param vars The variables object representing the current state of the
     *             numerical solution (primary and possibly secondary variables).
     */
    void solve(Variables& vars) override
    {
        const bool converged = solve_(vars);
        if (!converged)
            DUNE_THROW(NumericalProblem,
                Fmt::format("Newton solver didn't converge after {} iterations.\n", numSteps_));
    }

    /*!
     * \brief Run the Newton method to solve a non-linear system.
     *        The solver is responsible for all the strategic decisions.
     * \param vars The variables object representing the current state of the
     *             numerical solution (primary and possibly secondary variables).
     * \post If converged, the `Variables` will represent the solution. If convergence
     *       fails, they are in some intermediate, undefined state.
     */
    bool apply(Variables& vars) override
    {
        return solve_(vars);
    }

    /*!
     * \brief Called before the Newton method is applied to an
     *        non-linear system of equations.
     *
     * \param initVars The variables representing the initial solution
     */
    virtual void newtonBegin(Variables& initVars)
    {
        numSteps_ = 0;

        if constexpr (hasPriVarsSwitch<PriVarSwitchVariables>)
        {
            if constexpr (assemblerExportsVariables)
                priVarSwitchAdapter_->initialize(Backend::dofs(initVars), initVars);
            else // this assumes assembly with solution (i.e. Variables=SolutionVector)
                priVarSwitchAdapter_->initialize(initVars, this->assembler().gridVariables());
        }


        const auto& initSol = Backend::dofs(initVars);

        // write the initial residual if a convergence writer was set
        if (convergenceWriter_)
        {
            this->assembler().assembleResidual(initVars);

            // dummy vector, there is no delta before solving the linear system
            ResidualVector delta = LinearAlgebraNativeBackend::zeros(Backend::size(initSol));
            convergenceWriter_->write(initSol, delta, this->assembler().residual());
        }

        if (enablePartialReassembly_)
        {
            partialReassembler_->resetColors();
            resizeDistanceFromLastLinearization_(initSol, distanceFromLastLinearization_);
        }
    }

    /*!
     * \brief Returns true if another iteration should be done.
     *
     * \param varsCurrentIter The variables of the current Newton iteration
     * \param converged if the Newton method's convergence criterion was met in this step
     */
    virtual bool newtonProceed(const Variables &varsCurrentIter, bool converged)
    {
        if (numSteps_ < minSteps_)
            return true;
        else if (converged)
            return false; // we are below the desired tolerance
        else if (numSteps_ >= maxSteps_)
        {
            // We have exceeded the allowed number of steps. If the
            // maximum relative shift was reduced by a factor of at least 4,
            // we proceed even if we are above the maximum number of steps.
            if (enableShiftCriterion_)
                return shift_*4.0 < lastShift_;
            else
                return reduction_*4.0 < lastReduction_;
        }

        return true;
    }

    /*!
     * \brief Indicates the beginning of a Newton iteration.
     */
    virtual void newtonBeginStep(const Variables& vars)
    {
        lastShift_ = shift_;
        if (numSteps_ == 0)
        {
            lastReduction_ = 1.0;
        }
        else
        {
            lastReduction_ = reduction_;
        }
    }

    /*!
     * \brief Assemble the linear system of equations \f$\mathbf{A}x - b = 0\f$.
     *
     * \param vars The current iteration's variables
     */
    virtual void assembleLinearSystem(const Variables& vars)
    {
        assembleLinearSystem_(this->assembler(), vars);

        if (enablePartialReassembly_)
            partialReassembler_->report(comm_, endIterMsgStream_);
    }

    /*!
     * \brief Solve the linear system of equations \f$\mathbf{A}x - b = 0\f$.
     *
     * Throws Dumux::NumericalProblem if the linear solver didn't
     * converge.
     *
     * If the linear solver doesn't accept multitype matrices we copy the matrix
     * into a 1x1 block BCRS matrix for solving.
     *
     * \param deltaU The difference between the current and the next solution
     */
    void solveLinearSystem(ResidualVector& deltaU)
    {
        bool converged = false;

        try
        {
            if (numSteps_ == 0)
                initialResidual_ = this->linearSolver().norm(this->assembler().residual());

            // solve by calling the appropriate implementation depending on whether the linear solver
            // is capable of handling MultiType matrices or not
            converged = solveLinearSystem_(deltaU);
        }
        catch (const Dune::Exception &e)
        {
            if (verbosity_ >= 1)
                std::cout << "Newton: Caught exception from the linear solver: \"" << e.what() << "\"\n";

            converged = false;
        }

        // make sure all processes converged
        int convergedRemote = converged;
        if (comm_.size() > 1)
            convergedRemote = comm_.min(converged);

        if (!converged)
        {
            DUNE_THROW(NumericalProblem, "Linear solver did not converge");
            ++numLinearSolverBreakdowns_;
        }
        else if (!convergedRemote)
        {
            DUNE_THROW(NumericalProblem, "Linear solver did not converge on a remote process");
            ++numLinearSolverBreakdowns_; // we keep correct count for process 0
        }
    }

    /*!
     * \brief Update the current solution with a delta vector.
     *
     * The error estimates required for the newtonConverged() and
     * newtonProceed() methods should be updated inside this method.
     *
     * Different update strategies, such as line search and chopped
     * updates can be implemented. The default behavior is just to
     * subtract deltaU from uLastIter, i.e.
     * \f[ u^{k+1} = u^k - \Delta u^k \f]
     *
     * \param vars The variables after the current iteration
     * \param uLastIter The solution vector after the last iteration
     * \param deltaU The delta as calculated from solving the linear
     *               system of equations. This parameter also stores
     *               the updated solution.
     */
    void newtonUpdate(Variables& vars,
                      const SolutionVector& uLastIter,
                      const ResidualVector& deltaU)
    {
        if (enableShiftCriterion_ || enablePartialReassembly_)
            newtonUpdateShift_(uLastIter, deltaU);

        if (enablePartialReassembly_) {
            // Determine the threshold 'eps' that is used for the partial reassembly.
            // Every entity where the primary variables exhibit a relative shift
            // summed up since the last linearization above 'eps' will be colored
            // red yielding a reassembly.
            // The user can provide three parameters to influence the threshold:
            // 'minEps' by 'Newton.ReassemblyMinThreshold' (1e-1*shiftTolerance_ default)
            // 'maxEps' by 'Newton.ReassemblyMaxThreshold' (1e2*shiftTolerance_ default)
            // 'omega'  by 'Newton.ReassemblyShiftWeight'  (1e-3 default)
            // The threshold is calculated from the currently achieved maximum
            // relative shift according to the formula
            // eps = max( minEps, min(maxEps, omega*shift) ).
            // Increasing/decreasing 'minEps' leads to less/more reassembly if
            // 'omega*shift' is small, i.e., for the last Newton iterations.
            // Increasing/decreasing 'maxEps' leads to less/more reassembly if
            // 'omega*shift' is large, i.e., for the first Newton iterations.
            // Increasing/decreasing 'omega' leads to more/less first and last
            // iterations in this sense.
            using std::max;
            using std::min;
            auto reassemblyThreshold = max(reassemblyMinThreshold_,
                                           min(reassemblyMaxThreshold_,
                                               shift_*reassemblyShiftWeight_));

            updateDistanceFromLastLinearization_(uLastIter, deltaU);
            partialReassembler_->computeColors(this->assembler(),
                                               distanceFromLastLinearization_,
                                               reassemblyThreshold);

            // set the discrepancy of the red entities to zero
            for (unsigned int i = 0; i < distanceFromLastLinearization_.size(); i++)
                if (partialReassembler_->dofColor(i) == EntityColor::red)
                    distanceFromLastLinearization_[i] = 0;
        }

        if (useLineSearch_)
            lineSearchUpdate_(vars, uLastIter, deltaU);

        else if (useChop_)
            choppedUpdate_(vars, uLastIter, deltaU);

        else
        {
            auto uCurrentIter = uLastIter;
            Backend::axpy(-1.0, deltaU, uCurrentIter);
            solutionChanged_(vars, uCurrentIter);

            if (enableResidualCriterion_)
                computeResidualReduction_(vars);
        }
    }

    /*!
     * \brief Indicates that one Newton iteration was finished.
     *
     * \param vars The variables after the current Newton iteration
     * \param uLastIter The solution at the beginning of the current Newton iteration
     */
    virtual void newtonEndStep(Variables &vars,
                               const SolutionVector &uLastIter)
    {
        if constexpr (hasPriVarsSwitch<PriVarSwitchVariables>)
        {
            if constexpr (assemblerExportsVariables)
                priVarSwitchAdapter_->invoke(Backend::dofs(vars), vars);
            else // this assumes assembly with solution (i.e. Variables=SolutionVector)
                priVarSwitchAdapter_->invoke(vars, this->assembler().gridVariables());
        }

        ++numSteps_;

        if (verbosity_ >= 1)
        {
            if (enableDynamicOutput_)
                std::cout << '\r'; // move cursor to beginning of line

            const auto width = Fmt::formatted_size("{}", maxSteps_);
            std::cout << Fmt::format("Newton iteration {:{}} done", numSteps_, width);

            if (enableShiftCriterion_)
                std::cout << Fmt::format(", maximum relative shift = {:.4e}", shift_);
            if (enableResidualCriterion_ && enableAbsoluteResidualCriterion_)
                std::cout << Fmt::format(", residual = {:.4e}", residualNorm_);
            else if (enableResidualCriterion_)
                std::cout << Fmt::format(", residual reduction = {:.4e}", reduction_);

            std::cout << endIterMsgStream_.str() << "\n";
        }
        endIterMsgStream_.str("");

        // When the Newton iterations are done: ask the model to check whether it makes sense
        // TODO: how do we realize this? -> do this here in the Newton solver
        // model_().checkPlausibility();
    }

    /*!
     * \brief Called if the Newton method ended
     *        (not known yet if we failed or succeeded)
     */
    virtual void newtonEnd()  {}

    /*!
     * \brief Returns true if the error of the solution is below the
     *        tolerance.
     */
    virtual bool newtonConverged() const
    {
        // in case the model has a priVar switch and some some primary variables
        // actually switched their state in the last iteration, enforce another iteration
        if (priVarSwitchAdapter_->switched())
            return false;

        if (enableShiftCriterion_ && !enableResidualCriterion_)
        {
            return shift_ <= shiftTolerance_;
        }
        else if (!enableShiftCriterion_ && enableResidualCriterion_)
        {
            if(enableAbsoluteResidualCriterion_)
                return residualNorm_ <= residualTolerance_;
            else
                return reduction_ <= reductionTolerance_;
        }
        else if (satisfyResidualAndShiftCriterion_)
        {
            if(enableAbsoluteResidualCriterion_)
                return shift_ <= shiftTolerance_
                        && residualNorm_ <= residualTolerance_;
            else
                return shift_ <= shiftTolerance_
                        && reduction_ <= reductionTolerance_;
        }
        else if(enableShiftCriterion_ && enableResidualCriterion_)
        {
            if(enableAbsoluteResidualCriterion_)
                return shift_ <= shiftTolerance_
                        || residualNorm_ <= residualTolerance_;
            else
                return shift_ <= shiftTolerance_
                        || reduction_ <= reductionTolerance_;
        }
        else
        {
            return shift_ <= shiftTolerance_
                    || reduction_ <= reductionTolerance_
                    || residualNorm_ <= residualTolerance_;
        }

        return false;
    }

    /*!
     * \brief Called if the Newton method broke down.
     * This method is called _after_ newtonEnd()
     */
    virtual void newtonFail(Variables& u) {}

    /*!
     * \brief Called if the Newton method ended successfully
     * This method is called _after_ newtonEnd()
     */
    virtual void newtonSucceed()  {}

    /*!
     * \brief output statistics / report
     */
    void report(std::ostream& sout = std::cout) const
    {
        sout << '\n'
             << "Newton statistics\n"
             << "----------------------------------------------\n"
             << "-- Total Newton iterations:            " << totalWastedIter_ + totalSucceededIter_ << '\n'
             << "-- Total wasted Newton iterations:     " << totalWastedIter_ << '\n'
             << "-- Total succeeded Newton iterations:  " << totalSucceededIter_ << '\n'
             << "-- Average iterations per solve:       " << std::setprecision(3) << double(totalSucceededIter_) / double(numConverged_) << '\n'
             << "-- Number of linear solver breakdowns: " << numLinearSolverBreakdowns_ << '\n'
             << std::endl;
    }

    /*!
     * \brief reset the statistics
     */
    void resetReport()
    {
        totalWastedIter_ = 0;
        totalSucceededIter_ = 0;
        numConverged_ = 0;
        numLinearSolverBreakdowns_ = 0;
    }

    /*!
     * \brief Report the options and parameters this Newton is configured with
     */
    void reportParams(std::ostream& sout = std::cout) const
    {
        sout << "\nNewton solver configured with the following options and parameters:\n";
        // options
        if (useLineSearch_) sout << " -- Newton.UseLineSearch = true\n";
        if (useChop_) sout << " -- Newton.EnableChop = true\n";
        if (enablePartialReassembly_) sout << " -- Newton.EnablePartialReassembly = true\n";
        if (enableAbsoluteResidualCriterion_) sout << " -- Newton.EnableAbsoluteResidualCriterion = true\n";
        if (enableShiftCriterion_) sout << " -- Newton.EnableShiftCriterion = true (relative shift convergence criterion)\n";
        if (enableResidualCriterion_) sout << " -- Newton.EnableResidualCriterion = true\n";
        if (satisfyResidualAndShiftCriterion_) sout << " -- Newton.SatisfyResidualAndShiftCriterion = true\n";
        // parameters
        if (enableShiftCriterion_) sout << " -- Newton.MaxRelativeShift = " << shiftTolerance_ << '\n';
        if (enableAbsoluteResidualCriterion_) sout << " -- Newton.MaxAbsoluteResidual = " << residualTolerance_ << '\n';
        if (enableResidualCriterion_) sout << " -- Newton.ResidualReduction = " << reductionTolerance_ << '\n';
        sout << " -- Newton.MinSteps = " << minSteps_ << '\n';
        sout << " -- Newton.MaxSteps = " << maxSteps_ << '\n';
        sout << " -- Newton.TargetSteps = " << targetSteps_ << '\n';
        if (enablePartialReassembly_)
        {
            sout << " -- Newton.ReassemblyMinThreshold = " << reassemblyMinThreshold_ << '\n';
            sout << " -- Newton.ReassemblyMaxThreshold = " << reassemblyMaxThreshold_ << '\n';
            sout << " -- Newton.ReassemblyShiftWeight = " << reassemblyShiftWeight_ << '\n';
        }
        sout << " -- Newton.RetryTimeStepReductionFactor = " << retryTimeStepReductionFactor_ << '\n';
        sout << " -- Newton.MaxTimeStepDivisions = " << maxTimeStepDivisions_ << '\n';
        sout << std::endl;
    }

    /*!
     * \brief Suggest a new time-step size based on the old time-step
     *        size.
     *
     * The default behavior is to suggest the old time-step size
     * scaled by the ratio between the target iterations and the
     * iterations required to actually solve the last time-step.
     */
    Scalar suggestTimeStepSize(Scalar oldTimeStep) const
    {
        // be aggressive reducing the time-step size but
        // conservative when increasing it. the rationale is
        // that we want to avoid failing in the next Newton
        // iteration which would require another linearization
        // of the problem.
        if (numSteps_ > targetSteps_) {
            Scalar percent = Scalar(numSteps_ - targetSteps_)/targetSteps_;
            return oldTimeStep/(1.0 + percent);
        }

        Scalar percent = Scalar(targetSteps_ - numSteps_)/targetSteps_;
        return oldTimeStep*(1.0 + percent/1.2);
    }

    /*!
     * \brief Specifies the verbosity level
     */
    void setVerbosity(int val)
    { verbosity_ = val; }

    /*!
     * \brief Return the verbosity level
     */
    int verbosity() const
    { return verbosity_ ; }

    /*!
     * \brief Returns the parameter group
     */
    const std::string& paramGroup() const
    { return paramGroup_; }

    /*!
     * \brief Attach a convergence writer to write out intermediate results after each iteration
     */
    void attachConvergenceWriter(std::shared_ptr<ConvergenceWriter> convWriter)
    { convergenceWriter_ = convWriter; }

    /*!
     * \brief Detach the convergence writer to stop the output
     */
    void detachConvergenceWriter()
    { convergenceWriter_ = nullptr; }

    /*!
     * \brief Return the factor for reducing the time step after a Newton iteration has failed
     */
    Scalar retryTimeStepReductionFactor() const
    { return retryTimeStepReductionFactor_; }

    /*!
     * \brief Set the factor for reducing the time step after a Newton iteration has failed
     */
    void setRetryTimeStepReductionFactor(const Scalar factor)
    { retryTimeStepReductionFactor_ = factor; }

protected:

    /*!
     * \brief Update solution-dependent quantities like grid variables after the solution has changed.
     * \todo TODO: In case we stop support for old-style grid variables / assemblers at one point,
     *             this would become obsolete as only the update call to the backend would remain.
     */
    virtual void solutionChanged_(Variables& vars, const SolutionVector& uCurrentIter)
    {
        Backend::update(vars, uCurrentIter);

        if constexpr (!assemblerExportsVariables)
            this->assembler().updateGridVariables(Backend::dofs(vars));
    }

    void computeResidualReduction_(const Variables& vars)
    {
        // we assume that the assembler works on solution vectors
        // if it doesn't export the variables type
        if constexpr (!assemblerExportsVariables)
            this->assembler().assembleResidual(Backend::dofs(vars));
        else
            this->assembler().assembleResidual(vars);

        residualNorm_ = this->linearSolver().norm(this->assembler().residual());

        reduction_ = residualNorm_/initialResidual_;
    }

    bool enableResidualCriterion() const
    { return enableResidualCriterion_; }

    //! optimal number of iterations we want to achieve
    int targetSteps_;
    //! minimum number of iterations we do
    int minSteps_;
    //! maximum number of iterations we do before giving up
    int maxSteps_;
    //! actual number of steps done so far
    int numSteps_;

    // residual criterion variables
    Scalar reduction_;
    Scalar residualNorm_;
    Scalar lastReduction_;
    Scalar initialResidual_;

    // shift criterion variables
    Scalar shift_;
    Scalar lastShift_;

    //! message stream to be displayed at the end of iterations
    std::ostringstream endIterMsgStream_;


private:

    /*!
     * \brief Run the Newton method to solve a non-linear system.
     *        The solver is responsible for all the strategic decisions.
     */
    bool solve_(Variables& vars)
    {
        try
        {
            // newtonBegin may manipulate the solution
            newtonBegin(vars);

            // the given solution is the initial guess
            auto uLastIter = Backend::dofs(vars);
            ResidualVector deltaU = LinearAlgebraNativeBackend::zeros(Backend::size(Backend::dofs(vars)));
            Detail::Newton::assign(deltaU, Backend::dofs(vars));

            // setup timers
            Dune::Timer assembleTimer(false);
            Dune::Timer solveTimer(false);
            Dune::Timer updateTimer(false);

            // execute the method as long as the solver thinks
            // that we should do another iteration
            bool converged = false;
            while (newtonProceed(vars, converged))
            {
                // notify the solver that we're about to start
                // a new iteration
                newtonBeginStep(vars);

                // make the current solution to the old one
                if (numSteps_ > 0)
                    uLastIter = Backend::dofs(vars);

                if (verbosity_ >= 1 && enableDynamicOutput_)
                    std::cout << "Assemble: r(x^k) = dS/dt + div F - q;   M = grad r"
                              << std::flush;

                ///////////////
                // assemble
                ///////////////

                // linearize the problem at the current solution
                assembleTimer.start();
                assembleLinearSystem(vars);
                assembleTimer.stop();

                ///////////////
                // linear solve
                ///////////////

                // Clear the current line using an ansi escape
                // sequence.  for an explanation see
                // http://en.wikipedia.org/wiki/ANSI_escape_code
                const char clearRemainingLine[] = { 0x1b, '[', 'K', 0 };

                if (verbosity_ >= 1 && enableDynamicOutput_)
                    std::cout << "\rSolve: M deltax^k = r"
                              << clearRemainingLine << std::flush;

                // solve the resulting linear equation system
                solveTimer.start();

                // set the delta vector to zero before solving the linear system!
                deltaU = 0;

                solveLinearSystem(deltaU);
                solveTimer.stop();

                ///////////////
                // update
                ///////////////
                if (verbosity_ >= 1 && enableDynamicOutput_)
                    std::cout << "\rUpdate: x^(k+1) = x^k - deltax^k"
                              << clearRemainingLine << std::flush;

                updateTimer.start();
                // update the current solution (i.e. uOld) with the delta
                // (i.e. u). The result is stored in u
                newtonUpdate(vars, uLastIter, deltaU);
                updateTimer.stop();

                // tell the solver that we're done with this iteration
                newtonEndStep(vars, uLastIter);

                // if a convergence writer was specified compute residual and write output
                if (convergenceWriter_)
                {
                    this->assembler().assembleResidual(vars);
                    convergenceWriter_->write(Backend::dofs(vars), deltaU, this->assembler().residual());
                }

                // detect if the method has converged
                converged = newtonConverged();
            }

            // tell solver we are done
            newtonEnd();

            // reset state if Newton failed
            if (!newtonConverged())
            {
                totalWastedIter_ += numSteps_;
                newtonFail(vars);
                return false;
            }

            totalSucceededIter_ += numSteps_;
            numConverged_++;

            // tell solver we converged successfully
            newtonSucceed();

            if (verbosity_ >= 1) {
                const auto elapsedTot = assembleTimer.elapsed() + solveTimer.elapsed() + updateTimer.elapsed();
                std::cout << Fmt::format("Assemble/solve/update time: {:.2g}({:.2f}%)/{:.2g}({:.2f}%)/{:.2g}({:.2f}%)\n",
                                         assembleTimer.elapsed(), 100*assembleTimer.elapsed()/elapsedTot,
                                         solveTimer.elapsed(), 100*solveTimer.elapsed()/elapsedTot,
                                         updateTimer.elapsed(), 100*updateTimer.elapsed()/elapsedTot);
            }
            return true;

        }
        catch (const NumericalProblem &e)
        {
            if (verbosity_ >= 1)
                std::cout << "Newton: Caught exception: \"" << e.what() << "\"\n";

            totalWastedIter_ += numSteps_;

            newtonFail(vars);
            return false;
        }
    }

    //! assembleLinearSystem_ for assemblers that support partial reassembly
    template<class A>
    auto assembleLinearSystem_(const A& assembler, const Variables& vars)
    -> typename std::enable_if_t<decltype(isValid(Detail::Newton::supportsPartialReassembly())(assembler))::value, void>
    {
        this->assembler().assembleJacobianAndResidual(vars, partialReassembler_.get());
    }

    //! assembleLinearSystem_ for assemblers that don't support partial reassembly
    template<class A>
    auto assembleLinearSystem_(const A& assembler, const Variables& vars)
    -> typename std::enable_if_t<!decltype(isValid(Detail::Newton::supportsPartialReassembly())(assembler))::value, void>
    {
        this->assembler().assembleJacobianAndResidual(vars);
    }

    /*!
     * \brief Update the maximum relative shift of the solution compared to
     *        the previous iteration. Overload for "normal" solution vectors.
     *
     * \param uLastIter The current iterative solution
     * \param deltaU The difference between the current and the next solution
     */
    virtual void newtonUpdateShift_(const SolutionVector &uLastIter,
                                    const ResidualVector &deltaU)
    {
        auto uNew = uLastIter;
        Backend::axpy(-1.0, deltaU, uNew);
        shift_ = Detail::Newton::maxRelativeShift<Scalar>(uLastIter, uNew);

        if (comm_.size() > 1)
            shift_ = comm_.max(shift_);
    }

    virtual void lineSearchUpdate_(Variables &vars,
                                   const SolutionVector &uLastIter,
                                   const ResidualVector &deltaU)
    {
        Scalar lambda = 1.0;
        auto uCurrentIter = uLastIter;

        while (true)
        {
            Backend::axpy(-lambda, deltaU, uCurrentIter);
            solutionChanged_(vars, uCurrentIter);

            computeResidualReduction_(vars);

            if (reduction_ < lastReduction_ || lambda <= lineSearchMinRelaxationFactor_)
            {
                endIterMsgStream_ << Fmt::format(", residual reduction {:.4e}->{:.4e}@lambda={:.4f}", lastReduction_, reduction_, lambda);
                return;
            }

            // try with a smaller update and reset solution
            lambda *= 0.5;
            uCurrentIter = uLastIter;
        }
    }

    //! \note method must update the gridVariables, too!
    virtual void choppedUpdate_(Variables& vars,
                                const SolutionVector& uLastIter,
                                const ResidualVector& deltaU)
    {
        DUNE_THROW(Dune::NotImplemented,
                   "Chopped Newton update strategy not implemented.");
    }

    /*!
     * \brief Solve the linear system of equations \f$\mathbf{A}x - b = 0\f$.
     *
     * Throws Dumux::NumericalProblem if the linear solver didn't converge.
     */
    virtual bool solveLinearSystem_(ResidualVector& deltaU)
    {
        assert(this->checkSizesOfSubMatrices(this->assembler().jacobian()) && "Matrix blocks have wrong sizes!");

        return this->linearSolver().solve(
            this->assembler().jacobian(),
            deltaU,
            this->assembler().residual()
        );
    }

    //! initialize the parameters by reading from the parameter tree
    void initParams_(const std::string& group = "")
    {
        useLineSearch_ = getParamFromGroup<bool>(group, "Newton.UseLineSearch", false);
        lineSearchMinRelaxationFactor_ = getParamFromGroup<Scalar>(group, "Newton.LineSearchMinRelaxationFactor", 0.125);
        useChop_ = getParamFromGroup<bool>(group, "Newton.EnableChop", false);
        if(useLineSearch_ && useChop_)
            DUNE_THROW(Dune::InvalidStateException, "Use either linesearch OR chop!");

        enableAbsoluteResidualCriterion_ = getParamFromGroup<bool>(group, "Newton.EnableAbsoluteResidualCriterion", false);
        enableShiftCriterion_ = getParamFromGroup<bool>(group, "Newton.EnableShiftCriterion", true);
        enableResidualCriterion_ = getParamFromGroup<bool>(group, "Newton.EnableResidualCriterion", false) || enableAbsoluteResidualCriterion_;
        satisfyResidualAndShiftCriterion_ = getParamFromGroup<bool>(group, "Newton.SatisfyResidualAndShiftCriterion", false);
        enableDynamicOutput_ = getParamFromGroup<bool>(group, "Newton.EnableDynamicOutput", true);

        if (!enableShiftCriterion_ && !enableResidualCriterion_)
        {
            DUNE_THROW(Dune::NotImplemented,
                       "at least one of NewtonEnableShiftCriterion or "
                       << "NewtonEnableResidualCriterion has to be set to true");
        }

        setMaxRelativeShift(getParamFromGroup<Scalar>(group, "Newton.MaxRelativeShift", 1e-8));
        setMaxAbsoluteResidual(getParamFromGroup<Scalar>(group, "Newton.MaxAbsoluteResidual", 1e-5));
        setResidualReduction(getParamFromGroup<Scalar>(group, "Newton.ResidualReduction", 1e-5));
        setTargetSteps(getParamFromGroup<int>(group, "Newton.TargetSteps", 10));
        setMinSteps(getParamFromGroup<int>(group, "Newton.MinSteps", 2));
        setMaxSteps(getParamFromGroup<int>(group, "Newton.MaxSteps", 18));

        enablePartialReassembly_ = getParamFromGroup<bool>(group, "Newton.EnablePartialReassembly", false);
        reassemblyMinThreshold_ = getParamFromGroup<Scalar>(group, "Newton.ReassemblyMinThreshold", 1e-1*shiftTolerance_);
        reassemblyMaxThreshold_ = getParamFromGroup<Scalar>(group, "Newton.ReassemblyMaxThreshold", 1e2*shiftTolerance_);
        reassemblyShiftWeight_ = getParamFromGroup<Scalar>(group, "Newton.ReassemblyShiftWeight", 1e-3);

        maxTimeStepDivisions_ = getParamFromGroup<std::size_t>(group, "Newton.MaxTimeStepDivisions", 10);
        retryTimeStepReductionFactor_ = getParamFromGroup<Scalar>(group, "Newton.RetryTimeStepReductionFactor", 0.5);

        verbosity_ = comm_.rank() == 0 ? getParamFromGroup<int>(group, "Newton.Verbosity", 2) : 0;
        numSteps_ = 0;

        // output a parameter report
        if (verbosity_ >= 2)
            reportParams();
    }

    template<class SolA, class SolB>
    void updateDistanceFromLastLinearization_(const SolA& u, const SolB& uDelta)
    {
        if constexpr (Dune::IsNumber<SolA>::value)
        {
            auto nextPriVars = u;
            nextPriVars -= uDelta;

            // add the current relative shift for this degree of freedom
            auto shift = Detail::Newton::maxRelativeShift<Scalar>(u, nextPriVars);
            distanceFromLastLinearization_[0] += shift;
        }
        else
        {
            for (std::size_t i = 0; i < u.size(); ++i)
            {
                const auto& currentPriVars(u[i]);
                auto nextPriVars(currentPriVars);
                nextPriVars -= uDelta[i];

                // add the current relative shift for this degree of freedom
                auto shift = Detail::Newton::maxRelativeShift<Scalar>(currentPriVars, nextPriVars);
                distanceFromLastLinearization_[i] += shift;
            }
        }
    }

    template<class ...ArgsA, class...ArgsB>
    void updateDistanceFromLastLinearization_(const Dune::MultiTypeBlockVector<ArgsA...>& uLastIter,
                                              const Dune::MultiTypeBlockVector<ArgsB...>& deltaU)
    {
        DUNE_THROW(Dune::NotImplemented, "Reassembly for MultiTypeBlockVector");
    }

    template<class Sol>
    void resizeDistanceFromLastLinearization_(const Sol& u, std::vector<Scalar>& dist)
    {
        dist.assign(Backend::size(u), 0.0);
    }

    template<class ...Args>
    void resizeDistanceFromLastLinearization_(const Dune::MultiTypeBlockVector<Args...>& u,
                                              std::vector<Scalar>& dist)
    {
        DUNE_THROW(Dune::NotImplemented, "Reassembly for MultiTypeBlockVector");
    }

    //! The communication object
    Communication comm_;

    //! the verbosity level
    int verbosity_;

    Scalar shiftTolerance_;
    Scalar reductionTolerance_;
    Scalar residualTolerance_;

    // time step control
    std::size_t maxTimeStepDivisions_;
    Scalar retryTimeStepReductionFactor_;

    // further parameters
    bool useLineSearch_;
    Scalar lineSearchMinRelaxationFactor_;
    bool useChop_;
    bool enableAbsoluteResidualCriterion_;
    bool enableShiftCriterion_;
    bool enableResidualCriterion_;
    bool satisfyResidualAndShiftCriterion_;
    bool enableDynamicOutput_;

    //! the parameter group for getting parameters from the parameter tree
    std::string paramGroup_;

    // infrastructure for partial reassembly
    bool enablePartialReassembly_;
    std::unique_ptr<Reassembler> partialReassembler_;
    std::vector<Scalar> distanceFromLastLinearization_;
    Scalar reassemblyMinThreshold_;
    Scalar reassemblyMaxThreshold_;
    Scalar reassemblyShiftWeight_;

    // statistics for the optional report
    std::size_t totalWastedIter_ = 0; //! Newton steps in solves that didn't converge
    std::size_t totalSucceededIter_ = 0; //! Newton steps in solves that converged
    std::size_t numConverged_ = 0; //! total number of converged solves
    std::size_t numLinearSolverBreakdowns_ = 0; //! total number of linear solves that failed

    //! the class handling the primary variable switch
    std::unique_ptr<PrimaryVariableSwitchAdapter> priVarSwitchAdapter_;

    //! convergence writer
    std::shared_ptr<ConvergenceWriter> convergenceWriter_ = nullptr;
};

} // end namespace Dumux

#endif
