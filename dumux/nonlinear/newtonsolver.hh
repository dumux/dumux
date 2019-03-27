// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Nonlinear
 * \brief Reference implementation of the Newton solver.
 *
 * Usually this solver should be sufficient.
 */
#ifndef DUMUX_NEWTON_SOLVER_HH
#define DUMUX_NEWTON_SOLVER_HH

#include <cmath>
#include <memory>
#include <iostream>
#include <type_traits>

#include <dune/common/timer.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpicollectivecommunication.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/std/type_traits.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/multitypeblockvector.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/typetraits/vector.hh>
#include <dumux/common/typetraits/isvalid.hh>
#include <dumux/linear/linearsolveracceptsmultitypematrix.hh>
#include <dumux/linear/matrixconverter.hh>
#include <dumux/assembly/partialreassembler.hh>

#include "newtonconvergencewriter.hh"

#include <dune/common/deprecated.hh>

namespace Dumux {
namespace Detail {

//! helper struct detecting if an assembler supports partial reassembly
struct supportsPartialReassembly
{
    template<class Assembler>
    auto operator()(Assembler&& a)
    -> decltype(a.assembleJacobianAndResidual(std::declval<const typename Assembler::ResidualType&>(),
                                              std::declval<const PartialReassembler<Assembler>*>()))
    {}
};

//! helper aliases to extract a primary variable switch from the VolumeVariables (if defined, yields int otherwise)
template<class Assembler>
using DetectPVSwitch = typename Assembler::GridVariables::VolumeVariables::PrimaryVariableSwitch;

template<class Assembler>
using GetPVSwitch = Dune::Std::detected_or<int, DetectPVSwitch, Assembler>;

} // end namespace Detail

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
          class Comm = Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator> >
class NewtonSolver
{
    using Scalar = typename Assembler::Scalar;
    using JacobianMatrix = typename Assembler::JacobianMatrix;
    using SolutionVector = typename Assembler::ResidualType;
    using ConvergenceWriter = ConvergenceWriterInterface<SolutionVector>;

    using PrimaryVariableSwitch = typename Detail::GetPVSwitch<Assembler>::type;
    using HasPriVarsSwitch = typename Detail::GetPVSwitch<Assembler>::value_t; // std::true_type or std::false_type
    static constexpr bool hasPriVarsSwitch() { return HasPriVarsSwitch::value; };

public:

    using Communication = Comm;

    /*!
     * \brief The Constructor
     */
    NewtonSolver(std::shared_ptr<Assembler> assembler,
                 std::shared_ptr<LinearSolver> linearSolver,
                 const Communication& comm = Dune::MPIHelper::getCollectiveCommunication(),
                 const std::string& paramGroup = "")
    : endIterMsgStream_(std::ostringstream::out)
    , assembler_(assembler)
    , linearSolver_(linearSolver)
    , comm_(comm)
    , paramGroup_(paramGroup)
    {
        initParams_(paramGroup);

        // set the linear system (matrix & residual) in the assembler
        assembler_->setLinearSystem();

        // set a different default for the linear solver residual reduction
        // within the Newton the linear solver doesn't need to solve too exact
        linearSolver_->setResidualReduction(getParamFromGroup<Scalar>(paramGroup, "LinearSolver.ResidualReduction", 1e-6));

        // initialize the partial reassembler
        if (enablePartialReassembly_)
            partialReassembler_ = std::make_unique<Reassembler>(*assembler_);

        if (hasPriVarsSwitch())
        {
            const int priVarSwitchVerbosity = getParamFromGroup<int>(paramGroup, "PrimaryVariableSwitch.Verbosity", 1);
            priVarSwitch_ = std::make_unique<PrimaryVariableSwitch>(priVarSwitchVerbosity);
        }
    }

    virtual ~NewtonSolver() {}

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
     */
    template<class TimeLoop>
    DUNE_DEPRECATED_MSG("Use attachConvergenceWriter(convWriter) and solve(x, *timeLoop) instead")
    void solve(SolutionVector& uCurrentIter, TimeLoop& timeLoop,
               std::shared_ptr<ConvergenceWriter> convWriter)
    {
        if (!convergenceWriter_)
            attachConvergenceWriter(convWriter);

        solve(uCurrentIter, timeLoop);
    }
    /*!
     * \brief Run the Newton method to solve a non-linear system.
     *        Does time step control when the Newton fails to converge
     */
    template<class TimeLoop>
    void solve(SolutionVector& uCurrentIter, TimeLoop& timeLoop)
    {
        if (assembler_->isStationaryProblem())
            DUNE_THROW(Dune::InvalidStateException, "Using time step control with stationary problem makes no sense!");

        // try solving the non-linear system
        for (std::size_t i = 0; i <= maxTimeStepDivisions_; ++i)
        {
            // linearize & solve
            const bool converged = solve_(uCurrentIter);

            if (converged)
                return;

            else if (!converged && i < maxTimeStepDivisions_)
            {
                // set solution to previous solution
                uCurrentIter = assembler_->prevSol();

                // reset the grid variables to the previous solution
                assembler_->resetTimeStep(uCurrentIter);

                if (verbose_)
                    std::cout << "Newton solver did not converge with dt = "
                              << timeLoop.timeStepSize() << " seconds. Retrying with time step of "
                              << timeLoop.timeStepSize() * retryTimeStepReductionFactor_ << " seconds\n";

                // try again with dt = dt * retryTimeStepReductionFactor_
                timeLoop.setTimeStepSize(timeLoop.timeStepSize() * retryTimeStepReductionFactor_);
            }

            else
            {
                DUNE_THROW(NumericalProblem, "Newton solver didn't converge after "
                                             << maxTimeStepDivisions_ << " time-step divisions. dt="
                                             << timeLoop.timeStepSize() << '\n');
            }
        }
    }

    /*!
     * \brief Run the Newton method to solve a non-linear system.
     *        The solver is responsible for all the strategic decisions.
     */
    void solve(SolutionVector& uCurrentIter)
    {
        const bool converged = solve_(uCurrentIter);
        if (!converged)
            DUNE_THROW(NumericalProblem, "Newton solver didn't converge after "
                                         << numSteps_ << " iterations.\n");
    }

    /*!
     * \brief Called before the Newton method is applied to an
     *        non-linear system of equations.
     *
     * \param u The initial solution
     */
    virtual void newtonBegin(SolutionVector& u)
    {
        numSteps_ = 0;
        initPriVarSwitch_(u, HasPriVarsSwitch{});

        // write the initial residual if a convergence writer was set
        if (convergenceWriter_)
        {
            assembler_->assembleResidual(u);
            SolutionVector delta(u);
            delta = 0; // dummy vector, there is no delta before solving the linear system
            convergenceWriter_->write(u, delta, assembler_->residual());
        }
    }

    /*!
     * \brief Returns true if another iteration should be done.
     *
     * \param uCurrentIter The solution of the current Newton iteration
     * \param converged if the Newton method's convergence criterion was met in this step
     */
    virtual bool newtonProceed(const SolutionVector &uCurrentIter, bool converged)
    {
        if (numSteps_ < 2)
            return true; // we always do at least two iterations
        else if (converged) {
            return false; // we are below the desired tolerance
        }
        else if (numSteps_ >= maxSteps_) {
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
    virtual void newtonBeginStep(const SolutionVector& u)
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
     * \param uCurrentIter The current iteration's solution vector
     */
    virtual void assembleLinearSystem(const SolutionVector& uCurrentIter)
    {
        assembleLinearSystem_(*assembler_, uCurrentIter);

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
    void solveLinearSystem(SolutionVector& deltaU)
    {
        auto& b = assembler_->residual();

        try
        {
            if (numSteps_ == 0)
            {
                Scalar norm2 = b.two_norm2();
                if (comm_.size() > 1)
                    norm2 = comm_.sum(norm2);

                using std::sqrt;
                initialResidual_ = sqrt(norm2);
            }

            // solve by calling the appropriate implementation depending on whether the linear solver
            // is capable of handling MultiType matrices or not
            bool converged = solveLinearSystem_(deltaU);

            // make sure all processes converged
            int convergedRemote = converged;
            if (comm_.size() > 1)
                convergedRemote = comm_.min(converged);

            if (!converged) {
                DUNE_THROW(NumericalProblem,
                           "Linear solver did not converge");
            }
            else if (!convergedRemote) {
                DUNE_THROW(NumericalProblem,
                           "Linear solver did not converge on a remote process");
            }
        }
        catch (const Dune::Exception &e) {
            // make sure all processes converged
            int converged = 0;
            if (comm_.size() > 1)
                converged = comm_.min(converged);

            NumericalProblem p;
            p.message(e.what());
            throw p;
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
     * \param uCurrentIter The solution vector after the current iteration
     * \param uLastIter The solution vector after the last iteration
     * \param deltaU The delta as calculated from solving the linear
     *               system of equations. This parameter also stores
     *               the updated solution.
     */
    void newtonUpdate(SolutionVector &uCurrentIter,
                      const SolutionVector &uLastIter,
                      const SolutionVector &deltaU)
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
            partialReassembler_->computeColors(*assembler_,
                                               distanceFromLastLinearization_,
                                               reassemblyThreshold);

            // set the discrepancy of the red entities to zero
            for (unsigned int i = 0; i < distanceFromLastLinearization_.size(); i++)
                if (partialReassembler_->dofColor(i) == EntityColor::red)
                    distanceFromLastLinearization_[i] = 0;
        }

        if (useLineSearch_)
            lineSearchUpdate_(uCurrentIter, uLastIter, deltaU);

        else if (useChop_)
            choppedUpdate_(uCurrentIter, uLastIter, deltaU);

        else
        {
            uCurrentIter = uLastIter;
            uCurrentIter -= deltaU;

            if (enableResidualCriterion_)
                computeResidualReduction_(uCurrentIter);

            else
            {
                // If we get here, the convergence criterion does not require
                // additional residual evaluations. Thus, the grid variables have
                // not yet been updated to the new uCurrentIter.
                assembler_->updateGridVariables(uCurrentIter);
            }
        }
    }

    /*!
     * \brief Indicates that one Newton iteration was finished.
     *
     * \param uCurrentIter The solution after the current Newton iteration
     * \param uLastIter The solution at the beginning of the current Newton iteration
     */
    virtual void newtonEndStep(SolutionVector &uCurrentIter,
                               const SolutionVector &uLastIter)
    {
        invokePriVarSwitch_(uCurrentIter, HasPriVarsSwitch{});

        ++numSteps_;

        if (verbose_)
        {
            auto width = std::to_string(maxSteps_).size();
            std::cout << "\rNewton iteration " << std::setw(width) << numSteps_ << " done";

            auto formatFlags = std::cout.flags();
            auto prec = std::cout.precision();
            std::cout << std::scientific << std::setprecision(3);

            if (enableShiftCriterion_)
                std::cout << ", maximum relative shift = " << shift_;
            if (enableResidualCriterion_ && enableAbsoluteResidualCriterion_)
                std::cout << ", residual = " << residualNorm_;
            else if (enableResidualCriterion_)
                std::cout << ", residual reduction = " << reduction_;

            std::cout.flags(formatFlags);
            std::cout.precision(prec);

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
        if (priVarsSwitchedInLastIteration_)
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
    virtual void newtonFail(SolutionVector& u) {}

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
        if (verbose_)
            sout << '\n'
                 << "Newton statistics\n"
                 << "----------------------------------------------\n"
                 << "-- Total Newton iterations:           " << totalWastedIter_ + totalSucceededIter_ << '\n'
                 << "-- Total wasted Newton iterations:    " << totalWastedIter_ << '\n'
                 << "-- Total succeeded Newton iterations: " << totalSucceededIter_ << '\n'
                 << "-- Average iterations per solve:      " << std::setprecision(3) << double(totalSucceededIter_) / double(numConverged_) << '\n'
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
     * \brief Specifies if the Newton method ought to be chatty.
     */
    void setVerbose(bool val)
    { verbose_ = val; }

    /*!
     * \brief Returns true if the Newton method ought to be chatty.
     */
    bool verbose() const
    { return verbose_ ; }

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

protected:

    void computeResidualReduction_(const SolutionVector &uCurrentIter)
    {
        residualNorm_ = assembler_->residualNorm(uCurrentIter);
        reduction_ = residualNorm_;
        reduction_ /= initialResidual_;
    }

    bool enableResidualCriterion() const
    { return enableResidualCriterion_; }

    const LinearSolver& linearSolver() const
    { return *linearSolver_; }

    LinearSolver& linearSolver()
    { return *linearSolver_; }

    const Assembler& assembler() const
    { return *assembler_; }

    Assembler& assembler()
    { return *assembler_; }

    /*!
     * \brief Initialize the privar switch, noop if there is no priVarSwitch
     */
    void initPriVarSwitch_(SolutionVector&, std::false_type) {}

    /*!
     * \brief Initialize the privar switch
     */
    void initPriVarSwitch_(SolutionVector& sol, std::true_type)
    {
        priVarSwitch_->reset(sol.size());
        priVarsSwitchedInLastIteration_ = false;

        const auto& problem = assembler_->problem();
        const auto& fvGridGeometry = assembler_->fvGridGeometry();
        auto& gridVariables = assembler_->gridVariables();
        priVarSwitch_->updateBoundary(problem, fvGridGeometry, gridVariables, sol);
    }

    /*!
     * \brief Switch primary variables if necessary, noop if there is no priVarSwitch
     */
    void invokePriVarSwitch_(SolutionVector&, std::false_type) {}

    /*!
     * \brief Switch primary variables if necessary
     */
    void invokePriVarSwitch_(SolutionVector& uCurrentIter, std::true_type)
    {
        // update the variable switch (returns true if the pri vars at at least one dof were switched)
        // for disabled grid variable caching
        const auto& fvGridGeometry = assembler_->fvGridGeometry();
        const auto& problem = assembler_->problem();
        auto& gridVariables = assembler_->gridVariables();

        // invoke the primary variable switch
        priVarsSwitchedInLastIteration_ = priVarSwitch_->update(uCurrentIter, gridVariables,
                                                                problem, fvGridGeometry);

        if (priVarsSwitchedInLastIteration_)
        {
            for (const auto& element : elements(fvGridGeometry.gridView()))
            {
                // if the volume variables are cached globally, we need to update those where the primary variables have been switched
                priVarSwitch_->updateSwitchedVolVars(problem, element, fvGridGeometry, gridVariables, uCurrentIter);

                // if the flux variables are cached globally, we need to update those where the primary variables have been switched
                priVarSwitch_->updateSwitchedFluxVarsCache(problem, element, fvGridGeometry, gridVariables, uCurrentIter);
            }
        }
    }

    //! optimal number of iterations we want to achieve
    int targetSteps_;
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
    bool solve_(SolutionVector& uCurrentIter)
    {
        // the given solution is the initial guess
        SolutionVector uLastIter(uCurrentIter);
        SolutionVector deltaU(uCurrentIter);

        Dune::Timer assembleTimer(false);
        Dune::Timer solveTimer(false);
        Dune::Timer updateTimer(false);

        if (enablePartialReassembly_)
        {
            partialReassembler_->resetColors();
            resizeDistanceFromLastLinearization_(uCurrentIter, distanceFromLastLinearization_);
        }

        try
        {
            newtonBegin(uCurrentIter);

            // execute the method as long as the solver thinks
            // that we should do another iteration
            bool converged = false;
            while (newtonProceed(uCurrentIter, converged))
            {
                // notify the solver that we're about to start
                // a new timestep
                newtonBeginStep(uCurrentIter);

                // make the current solution to the old one
                if (numSteps_ > 0)
                    uLastIter = uCurrentIter;

                if (verbose_) {
                    std::cout << "Assemble: r(x^k) = dS/dt + div F - q;   M = grad r"
                              << std::flush;
                }

                ///////////////
                // assemble
                ///////////////

                // linearize the problem at the current solution
                assembleTimer.start();
                assembleLinearSystem(uCurrentIter);
                assembleTimer.stop();

                ///////////////
                // linear solve
                ///////////////

                // Clear the current line using an ansi escape
                // sequence.  for an explanation see
                // http://en.wikipedia.org/wiki/ANSI_escape_code
                const char clearRemainingLine[] = { 0x1b, '[', 'K', 0 };

                if (verbose_) {
                    std::cout << "\rSolve: M deltax^k = r";
                    std::cout << clearRemainingLine
                              << std::flush;
                }

                // solve the resulting linear equation system
                solveTimer.start();

                // set the delta vector to zero before solving the linear system!
                deltaU = 0;

                solveLinearSystem(deltaU);
                solveTimer.stop();

                ///////////////
                // update
                ///////////////
                if (verbose_) {
                    std::cout << "\rUpdate: x^(k+1) = x^k - deltax^k";
                    std::cout << clearRemainingLine;
                    std::cout.flush();
                }

                updateTimer.start();
                // update the current solution (i.e. uOld) with the delta
                // (i.e. u). The result is stored in u
                newtonUpdate(uCurrentIter, uLastIter, deltaU);
                updateTimer.stop();

                // tell the solver that we're done with this iteration
                newtonEndStep(uCurrentIter, uLastIter);

                // if a convergence writer was specified compute residual and write output
                if (convergenceWriter_)
                {
                    assembler_->assembleResidual(uCurrentIter);
                    convergenceWriter_->write(uCurrentIter, deltaU, assembler_->residual());
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
                newtonFail(uCurrentIter);
                return false;
            }

            totalSucceededIter_ += numSteps_;
            numConverged_++;

            // tell solver we converged successfully
            newtonSucceed();

            if (verbose_) {
                const auto elapsedTot = assembleTimer.elapsed() + solveTimer.elapsed() + updateTimer.elapsed();
                std::cout << "Assemble/solve/update time: "
                          <<  assembleTimer.elapsed() << "(" << 100*assembleTimer.elapsed()/elapsedTot << "%)/"
                          <<  solveTimer.elapsed() << "(" << 100*solveTimer.elapsed()/elapsedTot << "%)/"
                          <<  updateTimer.elapsed() << "(" << 100*updateTimer.elapsed()/elapsedTot << "%)"
                          << "\n";
            }
            return true;

        }
        catch (const NumericalProblem &e)
        {
            if (verbose_)
                std::cout << "Newton: Caught exception: \"" << e.what() << "\"\n";

            totalWastedIter_ += numSteps_;
            newtonFail(uCurrentIter);
            return false;
        }
    }

    //! assembleLinearSystem_ for assemblers that support partial reassembly
    template<class A>
    auto assembleLinearSystem_(const A& assembler, const SolutionVector& uCurrentIter)
    -> typename std::enable_if_t<decltype(isValid(Detail::supportsPartialReassembly())(assembler))::value, void>
    {
        assembler_->assembleJacobianAndResidual(uCurrentIter, partialReassembler_.get());
    }

    //! assembleLinearSystem_ for assemblers that don't support partial reassembly
    template<class A>
    auto assembleLinearSystem_(const A& assembler, const SolutionVector& uCurrentIter)
    -> typename std::enable_if_t<!decltype(isValid(Detail::supportsPartialReassembly())(assembler))::value, void>
    {
        assembler_->assembleJacobianAndResidual(uCurrentIter);
    }

    /*!
     * \brief Update the maximum relative shift of the solution compared to
     *        the previous iteration. Overload for "normal" solution vectors.
     *
     * \param uLastIter The current iterative solution
     * \param deltaU The difference between the current and the next solution
     */
    virtual void newtonUpdateShift_(const SolutionVector &uLastIter,
                                    const SolutionVector &deltaU)
    {
        shift_ = 0;
        newtonUpdateShiftImpl_(uLastIter, deltaU);

        if (comm_.size() > 1)
            shift_ = comm_.max(shift_);
    }

    template<class SolVec>
    void newtonUpdateShiftImpl_(const SolVec &uLastIter,
                                const SolVec &deltaU)
    {
        for (int i = 0; i < int(uLastIter.size()); ++i) {
            typename SolVec::block_type uNewI = uLastIter[i];
            uNewI -= deltaU[i];

            Scalar shiftAtDof = relativeShiftAtDof_(uLastIter[i], uNewI);
            using std::max;
            shift_ = max(shift_, shiftAtDof);
        }
    }

    template<class ...Args>
    void newtonUpdateShiftImpl_(const Dune::MultiTypeBlockVector<Args...> &uLastIter,
                                const Dune::MultiTypeBlockVector<Args...> &deltaU)
    {
        // There seems to be a bug in g++5 which which prevents compilation when
        // passing the call to the implementation directly to Dune::Hybrid::forEach.
        // We therefore store this call in a lambda and pass it to the for loop afterwards.
        auto doUpdate = [&](const auto subVectorIdx)
        {
            this->newtonUpdateShiftImpl_(uLastIter[subVectorIdx], deltaU[subVectorIdx]);
        };

        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(uLastIter)), doUpdate);
    }

    virtual void lineSearchUpdate_(SolutionVector &uCurrentIter,
                                   const SolutionVector &uLastIter,
                                   const SolutionVector &deltaU)
    {
        Scalar lambda = 1.0;
        SolutionVector tmp(uLastIter);

        while (true)
        {
            uCurrentIter = deltaU;
            uCurrentIter *= -lambda;
            uCurrentIter += uLastIter;

            computeResidualReduction_(uCurrentIter);

            if (reduction_ < lastReduction_ || lambda <= 0.125) {
                endIterMsgStream_ << ", residual reduction " << lastReduction_ << "->"  << reduction_ << "@lambda=" << lambda;
                return;
            }

            // try with a smaller update
            lambda /= 2.0;
        }
    }

    //! \note method must update the gridVariables, too!
    virtual void choppedUpdate_(SolutionVector &uCurrentIter,
                                const SolutionVector &uLastIter,
                                const SolutionVector &deltaU)
    {
        DUNE_THROW(Dune::NotImplemented,
                   "Chopped Newton update strategy not implemented.");
    }

    virtual bool solveLinearSystem_(SolutionVector& deltaU)
    {
        return solveLinearSystemImpl_(*linearSolver_,
                                      assembler_->jacobian(),
                                      deltaU,
                                      assembler_->residual());
    }

    /*!
     * \brief Solve the linear system of equations \f$\mathbf{A}x - b = 0\f$.
     *
     * Throws Dumux::NumericalProblem if the linear solver didn't
     * converge.
     *
     * Specialization for linear solvers that can handle MultiType matrices.
     *
     */
    template<class V = SolutionVector>
    typename std::enable_if_t<!isMultiTypeBlockVector<V>(), bool>
    solveLinearSystemImpl_(LinearSolver& ls,
                           JacobianMatrix& A,
                           SolutionVector& x,
                           SolutionVector& b)
    {
        //! Copy into a standard block vector.
        //! This is necessary for all model _not_ using a FieldVector<Scalar, blockSize> as
        //! primary variables vector in combination with UMFPack or SuperLU as their interfaces are hard coded
        //! to this field vector type in Dune ISTL
        //! Could be avoided for vectors that already have the right type using SFINAE
        //! but it shouldn't impact performance too much
        constexpr auto blockSize = JacobianMatrix::block_type::rows;
        using BlockType = Dune::FieldVector<Scalar, blockSize>;
        Dune::BlockVector<BlockType> xTmp; xTmp.resize(b.size());
        Dune::BlockVector<BlockType> bTmp(xTmp);
        for (unsigned int i = 0; i < b.size(); ++i)
            for (unsigned int j = 0; j < blockSize; ++j)
                bTmp[i][j] = b[i][j];

        const int converged = ls.solve(A, xTmp, bTmp);

        for (unsigned int i = 0; i < x.size(); ++i)
            for (unsigned int j = 0; j < blockSize; ++j)
                x[i][j] = xTmp[i][j];

        return converged;
    }


    /*!
     * \brief Solve the linear system of equations \f$\mathbf{A}x - b = 0\f$.
     *
     * Throws Dumux::NumericalProblem if the linear solver didn't
     * converge.
     *
     * Specialization for linear solvers that can handle MultiType matrices.
     *
     */

    template<class LS = LinearSolver, class V = SolutionVector>
    typename std::enable_if_t<linearSolverAcceptsMultiTypeMatrix<LS>() &&
                              isMultiTypeBlockVector<V>(), bool>
    solveLinearSystemImpl_(LinearSolver& ls,
                           JacobianMatrix& A,
                           SolutionVector& x,
                           SolutionVector& b)
    {
        // check matrix sizes
        assert(checkMatrix_(A) && "Sub blocks of MultiType matrix have wrong sizes!");

        // TODO: automatically derive the precondBlockLevel
        return ls.template solve</*precondBlockLevel=*/2>(A, x, b);
    }

    /*!
     * \brief Solve the linear system of equations \f$\mathbf{A}x - b = 0\f$.
     *
     * Throws Dumux::NumericalProblem if the linear solver didn't
     * converge.
     *
     * Specialization for linear solvers that cannot handle MultiType matrices.
     * We copy the matrix into a 1x1 block BCRS matrix before solving.
     *
     */
    template<class LS = LinearSolver, class V = SolutionVector>
    typename std::enable_if_t<!linearSolverAcceptsMultiTypeMatrix<LS>() &&
                              isMultiTypeBlockVector<V>(), bool>
    solveLinearSystemImpl_(LinearSolver& ls,
                           JacobianMatrix& A,
                           SolutionVector& x,
                           SolutionVector& b)
    {
        // check matrix sizes
        assert(checkMatrix_(A) && "Sub blocks of MultiType matrix have wrong sizes!");

        // create the bcrs matrix the IterativeSolver backend can handle
        const auto M = MatrixConverter<JacobianMatrix>::multiTypeToBCRSMatrix(A);

        // get the new matrix sizes
        const std::size_t numRows = M.N();
        assert(numRows == M.M());

        // create the vector the IterativeSolver backend can handle
        const auto bTmp = VectorConverter<SolutionVector>::multiTypeToBlockVector(b);
        assert(bTmp.size() == numRows);

        // create a blockvector to which the linear solver writes the solution
        using VectorBlock = typename Dune::FieldVector<Scalar, 1>;
        using BlockVector = typename Dune::BlockVector<VectorBlock>;
        BlockVector y(numRows);

        // solve
        const bool converged = ls.solve(M, y, bTmp);

        // copy back the result y into x
        if(converged)
            VectorConverter<SolutionVector>::retrieveValues(x, y);

        return converged;
    }

    //! helper method to assure the MultiType matrix's sub blocks have the correct sizes
    template<class M = JacobianMatrix>
    typename std::enable_if_t<!isBCRSMatrix<M>(), bool>
    checkMatrix_(const JacobianMatrix& A)
    {
        bool matrixHasCorrectSize = true;
        using namespace Dune::Hybrid;
        using namespace Dune::Indices;
        forEach(A, [&matrixHasCorrectSize](const auto& rowOfMultiTypeMatrix)
        {
            const auto numRowsLeftMostBlock = rowOfMultiTypeMatrix[_0].N();

            forEach(rowOfMultiTypeMatrix, [&matrixHasCorrectSize, &numRowsLeftMostBlock](const auto& subBlock)
            {
                if (subBlock.N() != numRowsLeftMostBlock)
                    matrixHasCorrectSize = false;
            });
        });
        return matrixHasCorrectSize;
    }

    //! initialize the parameters by reading from the parameter tree
    void initParams_(const std::string& group = "")
    {
        useLineSearch_ = getParamFromGroup<bool>(group, "Newton.UseLineSearch");
        useChop_ = getParamFromGroup<bool>(group, "Newton.EnableChop");
        if(useLineSearch_ && useChop_)
            DUNE_THROW(Dune::InvalidStateException, "Use either linesearch OR chop!");

        enableAbsoluteResidualCriterion_ = getParamFromGroup<bool>(group, "Newton.EnableAbsoluteResidualCriterion");
        enableShiftCriterion_ = getParamFromGroup<bool>(group, "Newton.EnableShiftCriterion");
        enableResidualCriterion_ = getParamFromGroup<bool>(group, "Newton.EnableResidualCriterion") || enableAbsoluteResidualCriterion_;
        satisfyResidualAndShiftCriterion_ = getParamFromGroup<bool>(group, "Newton.SatisfyResidualAndShiftCriterion");

        if (!enableShiftCriterion_ && !enableResidualCriterion_)
        {
            DUNE_THROW(Dune::NotImplemented,
                       "at least one of NewtonEnableShiftCriterion or "
                       << "NewtonEnableResidualCriterion has to be set to true");
        }

        setMaxRelativeShift(getParamFromGroup<Scalar>(group, "Newton.MaxRelativeShift"));
        setMaxAbsoluteResidual(getParamFromGroup<Scalar>(group, "Newton.MaxAbsoluteResidual"));
        setResidualReduction(getParamFromGroup<Scalar>(group, "Newton.ResidualReduction"));
        setTargetSteps(getParamFromGroup<int>(group, "Newton.TargetSteps"));
        setMaxSteps(getParamFromGroup<int>(group, "Newton.MaxSteps"));

        enablePartialReassembly_ = getParamFromGroup<bool>(group, "Newton.EnablePartialReassembly");
        reassemblyMinThreshold_ = getParamFromGroup<Scalar>(group, "Newton.ReassemblyMinThreshold", 1e-1*shiftTolerance_);
        reassemblyMaxThreshold_ = getParamFromGroup<Scalar>(group, "Newton.ReassemblyMaxThreshold", 1e2*shiftTolerance_);
        reassemblyShiftWeight_ = getParamFromGroup<Scalar>(group, "Newton.ReassemblyShiftWeight", 1e-3);

        maxTimeStepDivisions_ = getParamFromGroup<std::size_t>(group, "Newton.MaxTimeStepDivisions", 10);
        retryTimeStepReductionFactor_ = getParamFromGroup<Scalar>(group, "Newton.RetryTimeStepReductionFactor", 0.5);

        verbose_ = comm_.rank() == 0;
        numSteps_ = 0;
    }

    template<class Sol>
    void updateDistanceFromLastLinearization_(const Sol& u, const Sol& uDelta)
    {
        for (size_t i = 0; i < u.size(); ++i) {
            const auto& currentPriVars(u[i]);
            auto nextPriVars(currentPriVars);
            nextPriVars -= uDelta[i];

            // add the current relative shift for this degree of freedom
            auto shift = relativeShiftAtDof_(currentPriVars, nextPriVars);
            distanceFromLastLinearization_[i] += shift;
        }
    }

    template<class ...Args>
    void updateDistanceFromLastLinearization_(const Dune::MultiTypeBlockVector<Args...>& uLastIter,
                                              const Dune::MultiTypeBlockVector<Args...>& deltaU)
    {
        DUNE_THROW(Dune::NotImplemented, "Reassembly for MultiTypeBlockVector");
    }

    template<class Sol>
    void resizeDistanceFromLastLinearization_(const Sol& u, std::vector<Scalar>& dist)
    {
        dist.assign(u.size(), 0.0);
    }

    template<class ...Args>
    void resizeDistanceFromLastLinearization_(const Dune::MultiTypeBlockVector<Args...>& u,
                                              std::vector<Scalar>& dist)
    {
        DUNE_THROW(Dune::NotImplemented, "Reassembly for MultiTypeBlockVector");
    }

    /*!
     * \brief Returns the maximum relative shift between two vectors of
     *        primary variables.
     *
     * \param priVars1 The first vector of primary variables
     * \param priVars2 The second vector of primary variables
     */
    template<class PrimaryVariables>
    Scalar relativeShiftAtDof_(const PrimaryVariables &priVars1,
                               const PrimaryVariables &priVars2) const
    {
        Scalar result = 0.0;
        using std::abs;
        using std::max;
        // iterate over all primary variables
        for (int j = 0; j < PrimaryVariables::dimension; ++j) {
            Scalar eqErr = abs(priVars1[j] - priVars2[j]);
            eqErr /= max<Scalar>(1.0,abs(priVars1[j] + priVars2[j])/2);

            result = max(result, eqErr);
        }
        return result;
    }

    std::shared_ptr<Assembler> assembler_;
    std::shared_ptr<LinearSolver> linearSolver_;

    //! The communication object
    Communication comm_;

    //! switches on/off verbosity
    bool verbose_;

    Scalar shiftTolerance_;
    Scalar reductionTolerance_;
    Scalar residualTolerance_;

    // time step control
    std::size_t maxTimeStepDivisions_;
    Scalar retryTimeStepReductionFactor_;

    // further parameters
    bool useLineSearch_;
    bool useChop_;
    bool enableAbsoluteResidualCriterion_;
    bool enableShiftCriterion_;
    bool enableResidualCriterion_;
    bool satisfyResidualAndShiftCriterion_;

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

    //! the class handling the primary variable switch
    std::unique_ptr<PrimaryVariableSwitch> priVarSwitch_;
    //! if we switched primary variables in the last iteration
    bool priVarsSwitchedInLastIteration_ = false;

    //! convergence writer
    std::shared_ptr<ConvergenceWriter> convergenceWriter_ = nullptr;
};

} // end namespace Dumux

#endif
