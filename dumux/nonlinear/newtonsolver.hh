// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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

#include <dune/common/timer.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpicollectivecommunication.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/multitypeblockvector.hh>
#include <dune/istl/matrixmatrix.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/typetraits/vector.hh>
#include <dumux/common/typetraits/isvalid.hh>
#include <dumux/linear/linearsolveracceptsmultitypematrix.hh>
#include <dumux/linear/matrixconverter.hh>
#include <dumux/assembly/partialreassembler.hh>

#include <dumux/freeflow/navierstokes/indices.hh>

#include "newtonconvergencewriter.hh"

namespace Dumux {

//! helper struct detecting if an assembler supports partial reassembly
struct supportsPartialReassembly
{
    template<class Assembler>
    auto operator()(Assembler&& a)
    -> decltype(a.assembleCoefficientMatrixAndRHS(std::declval<const typename Assembler::ResidualType&>(),
                                              std::declval<const PartialReassembler<Assembler>*>()))
    {}
};

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
    using FVGridGeometry = typename Assembler::FVGridGeometry;
    using IndexType = typename FVGridGeometry::GridView::IndexSet::IndexType;
    using SolutionVector = typename Assembler::ResidualType;
    using ConvergenceWriter = ConvergenceWriterInterface<SolutionVector>;

    static constexpr auto faceIdx = FVGridGeometry::faceIdx();
    static constexpr auto cellCenterIdx = FVGridGeometry::cellCenterIdx();

    using CellCenterSolutionVector = typename Assembler::CellCenterSolutionVector;
    using FaceSolutionVector = typename Assembler::FaceSolutionVector;

    using SubControlVolume = typename FVGridGeometry::SubControlVolume;

    using CCToCCMatrixBlock = typename Assembler::CCToCCMatrixBlock;
    using FaceToFaceMatrixBlock = typename Assembler::FaceToFaceMatrixBlock;
    using FaceToCCMatrixBlock = typename Assembler::FaceToCCMatrixBlock;
    using CCToFaceMatrixBlock = typename Assembler::CCToFaceMatrixBlock;

    //TODO get pressureIdx from Indices instead (my idea was to have //     template<std::size_t id>
    //     using Indices = typename LocalResidual<id>::ModelTraits::Indices; in multidomainfvassembler and
    //     using Indices = typename ParentType::Indices<Dune::index_constant<0>()>; in staggeredfvassembler assembleCoefficientMatrixAndRHS
    //     enum {
    //         pressureIdx = Indices::pressureIdx
    //     };
    // as well as
    //     using Indices = typename Assembler::Indices;
    //here. However, the domainId is still a problem.
    static constexpr int pressureIdx = FVGridGeometry::GridView::dimension;

public:

    using Communication = Comm;

    /*!
     * \brief Constructor for stationary problems
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
    void solve(SolutionVector& uCurrentIter, TimeLoop& timeLoop,
               std::shared_ptr<ConvergenceWriter> convWriter = nullptr)
    {
        if (assembler_->isStationaryProblem())
            DUNE_THROW(Dune::InvalidStateException, "Using time step control with stationary problem makes no sense!");

        // try solving the non-linear system
        for (std::size_t i = 0; i <= maxTimeStepDivisions_; ++i)
        {
            // linearize & solve
            const bool converged = solve_(uCurrentIter, convWriter);

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
                              << timeLoop.timeStepSize()/2 << " seconds\n";

                // try again with dt = dt/2
                timeLoop.setTimeStepSize(timeLoop.timeStepSize()/2);
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
    void solve(SolutionVector& uCurrentIter, std::shared_ptr<ConvergenceWriter> convWriter = nullptr)
    {
        const bool converged = solve_(uCurrentIter, convWriter);
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
    virtual void newtonBegin(const SolutionVector& u)
    {
        numSteps_ = 0;
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
     * \param assembler The jacobian assembler
     * \param uCurrentIter The current iteration's solution vector
     */
    virtual void assembleLinearSystem(const SolutionVector& uCurrentIter)
    {
        assembleLinearSystem_(*assembler_, uCurrentIter);

        if (enablePartialReassembly_)
            partialReassembler_->report(comm_, endIterMsgStream_);
    }

    /*!
     * \brief Assemble the residual.
     *
     * \param assembler The assembler
     * \param uCurrentIter The current iteration's solution vector
     */
    virtual void assembleResidual(const SolutionVector& uCurrentIter)
    {
        assembleResidual_(*assembler_, uCurrentIter);
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
     * \param ls the linear solver
     * \param A The matrix of the linear system of equations
     * \param x The vector which solves the linear system
     * \param b The right hand side of the linear system
     */
    template<class Matrix, class Vector>
    void solveLinearSystem(const Matrix& A,
                           Vector& x,
                           const Vector& b)
    {
        //This is the residual that has been assembled when assembleResidual() has been called.
        //Hence, it is for uLastIter and not for the vector x or something like that.
        auto& res = assembler_->residual();

        try
        {
            if (numSteps_ == 0)
            {
                Scalar norm2 = res.two_norm2();
                if (comm_.size() > 1)
                    norm2 = comm_.sum(norm2);

                using std::sqrt;
                initialResidual_ = sqrt(norm2);
            }

            // solve by calling the appropriate implementation depending on whether the linear solver
            // is capable of handling MultiType matrices or not
            bool converged = solveLinearSystem_(A, x, b);

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
     * \param assembler The assembler (needed for global residual evaluation)
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
     * \param assembler The jacobian assembler
     * \param uCurrentIter The solution after the current Newton iteration
     * \param uLastIter The solution at the beginning of the current Newton iteration
     */
    virtual void newtonEndStep(SolutionVector &uCurrentIter,
                               const SolutionVector &uLastIter)
    {
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
    virtual void newtonEnd()  {
//         Scalar velocityUnderrelaxationFactor = getParamFromGroup<Scalar>("", "Underrelaxation.HundredTimesVelocityUnderrelaxationFactor")/100.;
//         Scalar pressureUnderrelaxationFactor = getParamFromGroup<Scalar>("", "Underrelaxation.HundredTimesPressureUnderrelaxationFactor")/100;
//
//         std::ofstream output("underrelaxation.txt", std::ios::app);
//
//         output << "vFactor = " << velocityUnderrelaxationFactor << ", pFactor = " << pressureUnderrelaxationFactor;
//
//         unsigned int algorithmType = getParamFromGroup<Scalar>("", "Algorithm.AlgorithmType", 0);
//         if (algorithmType == 1/*PISO*/){
//             Scalar secondStepVelocityUnderrelaxationFactor = getParamFromGroup<Scalar>("", "Underrelaxation.HundredTimesPISOSecondStepVelocityUnderrelaxationFactor")/100;
//             Scalar secondStepPressureUnderrelaxationFactor = getParamFromGroup<Scalar>("", "Underrelaxation.HundredTimesPISOSecondStepPressureUnderrelaxationFactor")/100;
//
//             output << ", PISOvFactor = " << secondStepVelocityUnderrelaxationFactor << ", PISOpFactor = " << secondStepPressureUnderrelaxationFactor;
//         }
//         output << ", StepsUntil 10^(-10) = " << numSteps_ << std::endl;
//
//         output.close();
    }

    /*!
     * \brief Returns true if the error of the solution is below the
     *        tolerance.
     */
    virtual bool newtonConverged() const
    {
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

protected:
    /*!
     * \brief Run the Newton method to solve a non-linear system.
     *        The solver is responsible for all the strategic decisions.
     */
    virtual bool solve_(SolutionVector& uCurrentIter, std::shared_ptr<ConvergenceWriter> convWriter = nullptr)
    {
        auto originalFullU = uCurrentIter;

        // make sure constructFullVectorFromReducedVector_ uses the correct boundary conditions even if the
        // initial conditions don't
        auto problem = (assembler_->problem());
        for (auto& scvIdx : problem.fixedPressureScvsIndexSet()){
            SubControlVolume scv = (assembler_-> fvGridGeometry()).scv(scvIdx);
            const auto dirichletAtCc = (assembler_->problem()).dirichletAtPos(scv.dofPosition());
            originalFullU[cellCenterIdx][scvIdx] = dirichletAtCc[pressureIdx];
        }
        //TODO add the same for originalFullU[faceIdx]. It is OK for the donea test, when starting with 0, because u=0 at the boundary, but it should be written in a more generic form!

        // the given solution is the initial guess
        SolutionVector uLastIter(uCurrentIter);

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
            while (newtonProceed(uCurrentIter, newtonConverged()))
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
                assembleResidual(uCurrentIter);
                assembleTimer.stop();

                // reduce uCurrentIter
                const auto boundaryScvfsIndexSet = (assembler_->problem()).dirichletBoundaryScvfsIndexSet();

                assembler_->removeSetOfEntriesFromVector(uCurrentIter[faceIdx], boundaryScvfsIndexSet);

                std::vector<IndexType> fixedPressureScvsIndexSet = (assembler_->problem()).fixedPressureScvsIndexSet();
                assembler_->removeSetOfEntriesFromVector(uCurrentIter[cellCenterIdx], fixedPressureScvsIndexSet);

//                 const auto& reducedULastIter = uCurrentIter;

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
                SolutionVector deltaU(uCurrentIter);
                deltaU = 0;

                solveLinearSystem(assembler_->reducedCoefficientMatrix(), deltaU, assembler_->reducedResidual());
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

                //full vectors from reduced ones
                //uCurrentIter

                uCurrentIter[faceIdx] = constructFullVectorFromReducedVector_(uCurrentIter[faceIdx], originalFullU[faceIdx], boundaryScvfsIndexSet);
                uCurrentIter[cellCenterIdx] = constructFullVectorFromReducedVector_(uCurrentIter[cellCenterIdx], originalFullU[cellCenterIdx], fixedPressureScvsIndexSet);

                // deltaU
                SolutionVector originalDeltaU;
                originalDeltaU[faceIdx].resize(originalFullU[faceIdx].size());
                originalDeltaU[cellCenterIdx].resize(originalFullU[cellCenterIdx].size());
                originalDeltaU = 0.;

                deltaU[faceIdx] = constructFullVectorFromReducedVector_(deltaU[faceIdx], originalDeltaU[faceIdx], boundaryScvfsIndexSet);
                deltaU[cellCenterIdx] = constructFullVectorFromReducedVector_(deltaU[cellCenterIdx], originalDeltaU[cellCenterIdx], fixedPressureScvsIndexSet);

                // update the current solution (i.e. uOld) with the delta
                // (i.e. u). The result is stored in u
                newtonUpdate(uCurrentIter, uLastIter, deltaU);
                updateTimer.stop();

                // tell the solver that we're done with this iteration
                newtonEndStep(uCurrentIter, uLastIter);

                // if a convergence writer was specified compute residual and write output
                if (convWriter)
                {
                    assembler_->assembleResidual(uCurrentIter);
                    convWriter->write(uLastIter, deltaU, assembler_->residual());
                }
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

    template<class VectorType, class IndexType>
    VectorType constructFullVectorFromReducedVector_(const VectorType& currentReducedVector,
                                                     const VectorType&  originalFullVector,
                                                     const std::vector<IndexType>& indices){
        if (indices.size() == 0) { return currentReducedVector; }

        std::vector<IndexType> tmpIndices = indices;
        std::sort (tmpIndices.begin(), tmpIndices.end());

        if (!(currentReducedVector.size() + indices.size() == originalFullVector.size())){
            std::cout << "Wrong sizes." << std::endl;
        }

        VectorType tmpVector;
        tmpVector.resize(originalFullVector.size());

        //fill intermediate reduced indices for A - delete rows
        int numBoundaryScvfsAlreadyHandled = 0;
        //k=0
        for (unsigned int i = 0; i < tmpIndices[0]; ++i){
            tmpVector[i] = currentReducedVector[i];
        }
        tmpVector[tmpIndices[0]] = originalFullVector[tmpIndices[0]];
        numBoundaryScvfsAlreadyHandled ++;
        for (unsigned int k = 1; k < tmpIndices.size(); ++k){
            //k is related to the boundary scvf up to which I want to go
            for (unsigned int i = (tmpIndices[k-1]+1); i < tmpIndices[k]; ++i){
                tmpVector[i]=currentReducedVector[i-numBoundaryScvfsAlreadyHandled];
            }
            numBoundaryScvfsAlreadyHandled ++;
            tmpVector[tmpIndices[k]] = originalFullVector[tmpIndices[k]];
        }
        for (unsigned int i = tmpIndices[tmpIndices.size()-1]+1; i < tmpVector.size(); ++i){
            tmpVector[i] = currentReducedVector[i-numBoundaryScvfsAlreadyHandled];
        }

        return tmpVector;
    }

    //! assembleLinearSystem_ for assemblers that support partial reassembly
    template<class A>
    auto assembleLinearSystem_(const A& assembler, const SolutionVector& uCurrentIter)
    -> typename std::enable_if_t<decltype(isValid(supportsPartialReassembly())(assembler))::value, void>
    {
        assembler_->assembleCoefficientMatrixAndRHS(uCurrentIter, partialReassembler_.get());
    }

    //! assembleLinearSystem_ for assemblers that don't support partial reassembly
    template<class A>
    auto assembleLinearSystem_(const A& assembler, const SolutionVector& uCurrentIter)
    -> typename std::enable_if_t<!decltype(isValid(supportsPartialReassembly())(assembler))::value, void>
    {
        assembler_->assembleCoefficientMatrixAndRHS(uCurrentIter);
    }

    template<class A>
    auto assembleResidual_(const A& assembler, const SolutionVector& uCurrentIter)
    {
        assembler_->assembleResidual(uCurrentIter);
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

    //underrelaxation
    virtual void lineSearchUpdate_(SolutionVector &uCurrentIter,
                                   const SolutionVector &uLastIter,
                                   const SolutionVector &deltaU)
    {
        //TODO adopt for SIMPLE
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

    template<class Vector, class Matrix>
    bool solveLinearSystem_(const Matrix& A,
                            Vector& x,
                            const Vector& b)
    {
        return solveLinearSystemImpl_(*linearSolver_,
                                      A,
                                      x,
                                      b);
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
    template<class Vector, class Matrix>
    typename std::enable_if_t<!isMultiTypeBlockVector<Vector>(), bool>
    solveLinearSystemImpl_(LinearSolver& ls,
                           const Matrix& A,
                           Vector& x,
                           const Vector& b)
    {
        //! Copy into a standard block vector.
        //! This is necessary for all model _not_ using a FieldVector<Scalar, blockSize> as
        //! primary variables vector in combination with UMFPack or SuperLU as their interfaces are hard coded
        //! to this field vector type in Dune ISTL
        //! Could be avoided for vectors that already have the right type using SFINAE
        //! but it shouldn't impact performance too much
        constexpr auto blockSize = Matrix::block_type::rows;
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

    template<class LS = LinearSolver, class Vector, class Matrix>
    typename std::enable_if_t<linearSolverAcceptsMultiTypeMatrix<LS>() &&
                              isMultiTypeBlockVector<Vector>(), bool>
    solveLinearSystemImpl_(LinearSolver& ls,
                           const Matrix& A,
                           Vector& x,
                           const Vector& b)
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
    template<class LS = LinearSolver, class Vector, class Matrix>
    typename std::enable_if_t<!linearSolverAcceptsMultiTypeMatrix<LS>() &&
                              isMultiTypeBlockVector<Vector>(), bool>
    solveLinearSystemImpl_(LinearSolver& ls,
                           const Matrix& A,
                           Vector& x,
                           const Vector& b)
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
    template<class Matrix>
    typename std::enable_if_t<!isBCRSMatrix<Matrix>(), bool>
    checkMatrix_(const Matrix& A)
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
};

template <class Assembler,
          class LinearSolver,
          class Reassembler = PartialReassembler<Assembler>,
          class Comm = Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator> >
class SimpleSolver : public NewtonSolver<Assembler, LinearSolver,
          Reassembler, Comm >
{
    using ParentType = NewtonSolver<Assembler, LinearSolver, Reassembler, Comm >;

    using Scalar = typename Assembler::Scalar;
    using FVGridGeometry = typename Assembler::FVGridGeometry;
    using IndexType = typename FVGridGeometry::GridView::IndexSet::IndexType;
    using SolutionVector = typename Assembler::ResidualType;
    using ConvergenceWriter = ConvergenceWriterInterface<SolutionVector>;

    static constexpr auto faceIdx = FVGridGeometry::faceIdx();
    static constexpr auto cellCenterIdx = FVGridGeometry::cellCenterIdx();

    using CellCenterSolutionVector = typename Assembler::CellCenterSolutionVector;
    using FaceSolutionVector = typename Assembler::FaceSolutionVector;

    using SubControlVolume = typename FVGridGeometry::SubControlVolume;

    using CCToCCMatrixBlock = typename Assembler::CCToCCMatrixBlock;
    using FaceToFaceMatrixBlock = typename Assembler::FaceToFaceMatrixBlock;
    using FaceToCCMatrixBlock = typename Assembler::FaceToCCMatrixBlock;
    using CCToFaceMatrixBlock = typename Assembler::CCToFaceMatrixBlock;

    static constexpr int pressureIdx = FVGridGeometry::GridView::dimension;

public:

    using Communication = Comm;

    /*!
     * \brief Constructor for stationary problems
     */
    SimpleSolver(std::shared_ptr<Assembler> assembler,
                 std::shared_ptr<LinearSolver> linearSolver,
                 const Communication& comm = Dune::MPIHelper::getCollectiveCommunication(),
                 const std::string& paramGroup = "")
    : ParentType (assembler, linearSolver, comm, paramGroup)
    {}

private:
    /*!
     * \brief Run the Newton method to solve a non-linear system.
     *        The solver is responsible for all the strategic decisions.
     */
    bool solve_(SolutionVector& uCurrentIter, std::shared_ptr<ConvergenceWriter> convWriter = nullptr)
    {
        auto originalFullU = uCurrentIter;

        // make sure constructFullVectorFromReducedVector_ uses the correct boundary conditions even if the
        // initial conditions don't
        auto problem = ((this->assembler_)->problem());
        for (auto& scvIdx : problem.fixedPressureScvsIndexSet()){
            SubControlVolume scv = ((this->assembler_)-> fvGridGeometry()).scv(scvIdx);
            const auto dirichletAtCc = ((this->assembler_)->problem()).dirichletAtPos(scv.dofPosition());
            originalFullU[cellCenterIdx][scvIdx] = dirichletAtCc[pressureIdx];
        }
        //TODO add the same for originalFullU[faceIdx]. It is OK for the donea test, when starting with 0, because u=0 at the boundary, but it should be written in a more generic form!

        // the given solution is the initial guess
        SolutionVector uLastIter(uCurrentIter);

        Dune::Timer assembleTimer(false);
        Dune::Timer solveTimer(false);
        Dune::Timer updateTimer(false);

        if (this->enablePartialReassembly_)
        {
            (this->partialReassembler_)->resetColors();
            this->resizeDistanceFromLastLinearization_(uCurrentIter, this->distanceFromLastLinearization_);
        }

        try
        {
            this->newtonBegin(uCurrentIter);

            // execute the method as long as the solver thinks
            // that we should do another iteration
            while (this->newtonProceed(uCurrentIter, this->newtonConverged()))
            {
                // notify the solver that we're about to start
                // a new timestep
                this->newtonBeginStep(uCurrentIter);

                // make the current solution to the old one
                if (this->numSteps_ > 0)
                    uLastIter = uCurrentIter;

                if (this->verbose_) {
                    std::cout << "Assemble: r(x^k) = dS/dt + div F - q;   M = grad r"
                              << std::flush;
                }

                ///////////////
                // assemble
                ///////////////

                // linearize the problem at the current solution
                assembleTimer.start();
                this->assembleLinearSystem(uCurrentIter);
                this->assembleResidual(uCurrentIter);
                assembleTimer.stop();

                // reduce uCurrentIter
                const auto boundaryScvfsIndexSet = ((this->assembler_)->problem()).dirichletBoundaryScvfsIndexSet();


                (this->assembler_)->removeSetOfEntriesFromVector(uCurrentIter[faceIdx], boundaryScvfsIndexSet);

                std::vector<IndexType> fixedPressureScvsIndexSet = ((this->assembler_)->problem()).fixedPressureScvsIndexSet();
                (this->assembler_)->removeSetOfEntriesFromVector(uCurrentIter[cellCenterIdx], fixedPressureScvsIndexSet);

                const auto& reducedULastIter = uCurrentIter;

                ///////////////
                // linear solve
                ///////////////

                // Clear the current line using an ansi escape
                // sequence.  for an explanation see
                // http://en.wikipedia.org/wiki/ANSI_escape_code
                const char clearRemainingLine[] = { 0x1b, '[', 'K', 0 };

                if (this->verbose_) {
                    std::cout << "\rSolve: M deltax^k = r";
                    std::cout << clearRemainingLine
                              << std::flush;
                }

                // solve the resulting linear equation system
                solveTimer.start();

                //To see which order of [faceIdx] and [cellCenterIdx] is correct, refer to staggeredlocalassembler.hh, function assembleCoefficientMatrixAndRHS, where the coefficientMatrix is filled.
                auto& A = (this->assembler_)->reducedCoefficientMatrix()[faceIdx][faceIdx];
                auto& B = (this->assembler_)->reducedCoefficientMatrix()[faceIdx][cellCenterIdx];
                auto& C = (this->assembler_)->reducedCoefficientMatrix()[cellCenterIdx][faceIdx];
                auto& rN = (this->assembler_)->reducedRHS()[faceIdx];
                auto& rC = (this->assembler_)->reducedRHS()[cellCenterIdx];

                const std::size_t numDofsFaceReduced = A.N();
                const std::size_t numDofsCellCenterReduced = B.M();
//                 Dune::writeVectorToMatlab(rC, "rC.m");
//                 Dune::writeVectorToMatlab(rN, "rN.m");
//                 Dune::writeMatrixToMatlab(A, "A.m");
//                 Dune::writeMatrixToMatlab(B, "B.m");
//                 Dune::writeMatrixToMatlab(C, "C.m");
//                 exit(0);

                FaceToCCMatrixBlock invAB;
                FaceToFaceMatrixBlock invDiagA;

                unsigned int algorithmType = getParamFromGroup<Scalar>("Algorithm", "Algorithm.AlgorithmType", 0);

                // [ Algorithm ]
                // # 0: SIMPLE
                // # 1: PISO
                // # 2: SIMPLER
                // # 3: SIMPLEC
                // # do not use SIMPLEC for Stokes
                // AlgorithmType = 2

                if(algorithmType == 0){
                //SIMPLE
                    invAB.setBuildMode(FaceToCCMatrixBlock::BuildMode::random);

                    Dune::MatrixIndexSet pattern(B.N(), B.M());

                    for (int i = 0; i < B.N(); ++i)
                    {
                        for (int j = 0; j < B.M(); ++j)
                        {
                            pattern.add(i,j);
                        }
                    }
                    pattern.exportIdx(invAB);

                    std::vector<Dune::BlockVector<Dune::FieldVector<Scalar,1>>> columnsOfB;
                    //B.M() is number of columns

                    columnsOfB.resize(B.M());

                    for (int i = 0; i < columnsOfB.size(); ++i)
                    {
                        columnsOfB[i].resize(B.N());
                        columnsOfB[i] = 0.0;
                    }

                    for (typename FaceToCCMatrixBlock::RowIterator i = B.begin(); i != B.end(); ++i)
                    {
                        for (typename FaceToCCMatrixBlock::ColIterator j = B[i.index()].begin(); j != B[i.index()].end(); ++j)
                        {
                            (columnsOfB[j.index()])[i.index()][0] = B[i.index()][j.index()][0][0];
                        }
                    }

    //                 for (const auto& vec : columnsOfB)
    //                 {
    //                     std::cout << "new columnvector ";
    //                     for (const auto& elem : vec)
    //                     {
    //                         std::cout << elem << ", ";
    //                     }
    //                     std::cout << std::endl;
    //                 }

                    for (int i = 0; i < B.M(); ++i)
                    {
                        Dune::BlockVector<Dune::FieldVector<Scalar,1>> ithColumnOfInvAB;
                        ithColumnOfInvAB.resize(B.N());
                        this->solveLinearSystem(A, ithColumnOfInvAB, columnsOfB[i]);
                        for (int j = 0; j < B.N(); ++j)
                        {
                            //[0] is the block size in the blockmatrix, compare in solveLinearSystemImpl_
                            invAB[j][i] = ithColumnOfInvAB[j][0];
                        }
                    }
                }
                else
                {
                    //SIMPLER, SIMPLEC or PISO
                    //get a diagonal matrix
                    //inspired by dumux/linear/amgbackend.hh
                    invDiagA.setBuildMode(FaceToFaceMatrixBlock::random);
                    setInvDiagAPattern_(invDiagA, numDofsFaceReduced);
                    typename FaceToFaceMatrixBlock::RowIterator row = A.begin();

                    for(; row != A.end(); ++row)
                    {
                        using size_type = typename FaceToFaceMatrixBlock::size_type;
                        size_type rowIdx = row.index();

                        //invDigaA = inverse(diagonal(A))
                        if (algorithmType == 3){
                            //SIMPLEC
                            Scalar rowSum = 0.0;
                            typename FaceToFaceMatrixBlock::ColIterator col = A[rowIdx].begin();
                            for (; col != A[rowIdx].end(); ++col)
                            {
                                size_type colIdx = col.index();
                                rowSum += A[rowIdx][colIdx];
                            }

                            invDiagA[rowIdx][rowIdx] = 1./rowSum;
                        }
                        else
                        {
                            invDiagA[rowIdx][rowIdx] = 1./(A[rowIdx][rowIdx]);
                        }
                    }

                    Dune::matMultMat(invAB, invDiagA, B);
                }

                CCToCCMatrixBlock matrixPressureStep;
                Dune::matMultMat(matrixPressureStep, C, invAB);

                if(algorithmType == 2){
                //SIMPLER
                    FaceToFaceMatrixBlock invDiagAAminusOne = getInvDiagAAminusOne_(invDiagA, A, numDofsFaceReduced);

                    FaceSolutionVector uHat;
                    uHat.resize(numDofsFaceReduced);
                    invDiagA.mv(rN, uHat);
                    invDiagAAminusOne.mmv(uCurrentIter[faceIdx], uHat);

                    auto resSIMPLERStep = rC;
                    resSIMPLERStep *= -1.;
                    C.umv(uHat, resSIMPLERStep);

                    this->solveLinearSystem(matrixPressureStep, uCurrentIter[cellCenterIdx], resSIMPLERStep);

                    //underrelaxation of pressure
                    Scalar pressureUnderrelaxationFactor = getParamFromGroup<Scalar>("Underrelaxation", "Underrelaxation.HundredTimesPressureUnderrelaxationFactor", 100)/100.;

                    //if anywhere, pressure has to be underrelaxed here. The reason is, that this pressure is only used in steps 3-5, and not in the next iteration step
                    auto oldP = reducedULastIter[cellCenterIdx];
                    oldP *= (1 - pressureUnderrelaxationFactor);
                    auto newP = uCurrentIter[cellCenterIdx];
                    newP *= pressureUnderrelaxationFactor;

                    oldP += newP;

                    uCurrentIter[cellCenterIdx] = oldP;
                }

                //velocity step
                auto velocityStepRHS = rN;
                velocityStepRHS *= -1.0;

                A.umv(uCurrentIter[faceIdx], velocityStepRHS);
                B.umv(uCurrentIter[cellCenterIdx], velocityStepRHS);

                FaceSolutionVector deltaUTilde;
                deltaUTilde.resize(velocityStepRHS.size());
                //TODO required (set zero before solving the linear system)!?
                deltaUTilde = 0;

                this->solveLinearSystem(A, deltaUTilde, velocityStepRHS);

                //pressure step
                auto pressureStepRHS = rC;
                C.mmv(uCurrentIter[faceIdx], pressureStepRHS);
                C.umv(deltaUTilde, pressureStepRHS);

                CellCenterSolutionVector pressureCorrection;
                pressureCorrection.resize(pressureStepRHS.size());
                //TODO required (set zero before solving the linear system)!?
                pressureCorrection = 0;

                this->solveLinearSystem(matrixPressureStep, pressureCorrection, pressureStepRHS);

                //calculate deltaU
                SolutionVector deltaU;
                //deltaU[faceIdx] = invAB * pressureCorrection
                deltaU[faceIdx] = deltaUTilde;
                invAB.mmv(pressureCorrection, deltaU[faceIdx]);
                deltaU[cellCenterIdx] = pressureCorrection;

                if (algorithmType != 2) {
                    Scalar pressureUnderrelaxationFactor = getParamFromGroup<Scalar>("Underrelaxation", "Underrelaxation.HundredTimesPressureUnderrelaxationFactor", 100)/100.;
                    deltaU[cellCenterIdx] *= pressureUnderrelaxationFactor;
                }

                Scalar velocityUnderrelaxationFactor = getParamFromGroup<Scalar>("Underrelaxation", "Underrelaxation.HundredTimesVelocityUnderrelaxationFactor", 100)/100.;
                deltaU[faceIdx] *= velocityUnderrelaxationFactor;

                solveTimer.stop();

                ///////////////
                // update
                ///////////////
                if (this->verbose_) {
                    std::cout << "\rUpdate: x^(k+1) = x^k - deltax^k";
                    std::cout << clearRemainingLine;
                    std::cout.flush();
                }

                updateTimer.start();

                //full vectors from reduced ones
                //uCurrentIter
                uCurrentIter[faceIdx] = this->constructFullVectorFromReducedVector_(uCurrentIter[faceIdx], originalFullU[faceIdx], boundaryScvfsIndexSet);
                uCurrentIter[cellCenterIdx] = this->constructFullVectorFromReducedVector_(uCurrentIter[cellCenterIdx], originalFullU[cellCenterIdx], fixedPressureScvsIndexSet);

                // deltaU
                SolutionVector originalDeltaU;
                originalDeltaU[faceIdx].resize(originalFullU[faceIdx].size());
                originalDeltaU[cellCenterIdx].resize(originalFullU[cellCenterIdx].size());
                originalDeltaU = 0.;
                deltaU[faceIdx] = this->constructFullVectorFromReducedVector_(deltaU[faceIdx], originalDeltaU[faceIdx], boundaryScvfsIndexSet);
                deltaU[cellCenterIdx] = this->constructFullVectorFromReducedVector_(deltaU[cellCenterIdx], originalDeltaU[cellCenterIdx], fixedPressureScvsIndexSet);

                if(algorithmType == 2){
                    //SIMPLER
                    deltaU[cellCenterIdx] = 0.;
                }

                // update the current solution (i.e. uOld) with the delta
                // (i.e. u). The result is stored in u
                this->newtonUpdate(uCurrentIter, uLastIter, deltaU);
                updateTimer.stop();

                if (algorithmType == 1){
                    //PISO
                    SolutionVector secondULastIter(uCurrentIter);

                    FaceToFaceMatrixBlock invDiagAAminusOne = getInvDiagAAminusOne_(invDiagA, A, numDofsFaceReduced);

                    CCToFaceMatrixBlock matrixForSecondPressureStepRHS;
                    Dune::matMultMat(matrixForSecondPressureStepRHS, C, invDiagAAminusOne);

                    // reduce uCurrentIter
                    (this->assembler_)->removeSetOfEntriesFromVector(uCurrentIter[faceIdx], boundaryScvfsIndexSet);
                    (this->assembler_)->removeSetOfEntriesFromVector(uCurrentIter[cellCenterIdx], fixedPressureScvsIndexSet);

                    // reduce uLastIter[faceIdx]
                    FaceSolutionVector uLastIterFaceReduced = uLastIter[faceIdx];
                    (this->assembler_)->removeSetOfEntriesFromVector(uLastIterFaceReduced, boundaryScvfsIndexSet);

                    CellCenterSolutionVector secondPressureStepRHS;
                    secondPressureStepRHS.resize(numDofsCellCenterReduced);
                    FaceSolutionVector deltaUForSecondPressureStepRHS = uCurrentIter[faceIdx];
                    deltaUForSecondPressureStepRHS -= uLastIterFaceReduced;
                    matrixForSecondPressureStepRHS.umv(deltaUForSecondPressureStepRHS, secondPressureStepRHS);

                    const auto& matrixSecondPressureStep = matrixPressureStep;

                    CellCenterSolutionVector secondPressureCorrection;
                    secondPressureCorrection.resize(secondPressureStepRHS.size());
                    //TODO required (set zero before solving the linear system)!?
                    secondPressureCorrection = 0;

                    this->solveLinearSystem(matrixSecondPressureStep, secondPressureCorrection, secondPressureStepRHS);

                    SolutionVector secondDeltaU;
                    secondDeltaU[faceIdx] = uLastIterFaceReduced;
                    secondDeltaU[faceIdx] -= uCurrentIter[faceIdx];
                    invDiagAAminusOne.mv(secondDeltaU[faceIdx], secondDeltaU[faceIdx]);

                    invAB.umv(secondPressureCorrection, secondDeltaU[faceIdx]);

                    secondDeltaU[cellCenterIdx] = secondPressureCorrection;

                    //underrelaxation
                    Scalar secondStepVelocityUnderrelaxationFactor = getParamFromGroup<Scalar>("Underrelaxation", "Underrelaxation.HundredTimesPISOSecondStepVelocityUnderrelaxationFactor", 100)/100.;
                    Scalar secondStepPressureUnderrelaxationFactor = getParamFromGroup<Scalar>("Underrelaxation", "Underrelaxation.HundredTimesPISOSecondStepPressureUnderrelaxationFactor", 100)/100.;

                    secondDeltaU[cellCenterIdx] *= secondStepPressureUnderrelaxationFactor;
                    secondDeltaU[faceIdx] *= secondStepVelocityUnderrelaxationFactor;

                    //full vectors from reduced ones
                    //uCurrentIter
                    uCurrentIter[faceIdx] = this->constructFullVectorFromReducedVector_(uCurrentIter[faceIdx], originalFullU[faceIdx], boundaryScvfsIndexSet);
                    uCurrentIter[cellCenterIdx] = this->constructFullVectorFromReducedVector_(uCurrentIter[cellCenterIdx], originalFullU[cellCenterIdx], fixedPressureScvsIndexSet);

                    // deltaU
                    secondDeltaU[faceIdx] = this->constructFullVectorFromReducedVector_(secondDeltaU[faceIdx], originalDeltaU[faceIdx], boundaryScvfsIndexSet);
                    secondDeltaU[cellCenterIdx] = this->constructFullVectorFromReducedVector_(secondDeltaU[cellCenterIdx], originalDeltaU[cellCenterIdx], fixedPressureScvsIndexSet);


                    this->newtonUpdate(uCurrentIter, secondULastIter, secondDeltaU);
                }

                // tell the solver that we're done with this iteration
                this->newtonEndStep(uCurrentIter, uLastIter);

                // if a convergence writer was specified compute residual and write output
                if (convWriter)
                {
                    (this->assembler_)->assembleResidual(uCurrentIter);
                    convWriter->write(uLastIter, deltaU, (this->assembler_)->residual());
                }
            }

            // tell solver we are done
            this->newtonEnd();

            // reset state if Newton failed
            if (!(this->newtonConverged()))
            {
                this->totalWastedIter_ += this->numSteps_;
                this->newtonFail(uCurrentIter);

                return false;
            }

            this->totalSucceededIter_ += this->numSteps_;
            this->numConverged_++;

            // tell solver we converged successfully
            this->newtonSucceed();

            if (this->verbose_) {
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
            if (this->verbose_)
                std::cout << "Newton: Caught exception: \"" << e.what() << "\"\n";

            this->totalWastedIter_ += this->numSteps_;
            this->newtonFail(uCurrentIter);
            return false;
        }
    }

    FaceToFaceMatrixBlock getInvDiagAAminusOne_(FaceToFaceMatrixBlock& A, FaceToFaceMatrixBlock& invDiagA, std::size_t numDofsFaceReduced){
        //get a unity matrix of size numDofsFaceReduced
        FaceToFaceMatrixBlock ones;
        ones.setBuildMode(FaceToFaceMatrixBlock::random);
        setOnesPattern_(ones, numDofsFaceReduced);
        typename FaceToFaceMatrixBlock::RowIterator row = A.begin();
        row = A.begin();
        for(; row != A.end(); ++row)
        {
            using size_type = typename FaceToFaceMatrixBlock::size_type;
            size_type rowIdx = row.index();

            ones[rowIdx][rowIdx] = 1.;
        }

        FaceToFaceMatrixBlock invDiagAAminusOne;
        Dune::matMultMat(invDiagAAminusOne, invDiagA, A);
        invDiagAAminusOne -= ones;

        return invDiagAAminusOne;
    }

    /*!
     * \brief Resizes the  matrix invDiagA and sets the matrix' sparsity pattern.
     */
    void setInvDiagAPattern_(FaceToFaceMatrixBlock& invDiagA, std::size_t numDofsFaceReduced)
    {
        // set the size of the sub-matrizes
        invDiagA.setSize(numDofsFaceReduced, numDofsFaceReduced);

        // set occupation pattern of the coefficient matrix
        Dune::MatrixIndexSet occupationPatternInvDiagA;
        occupationPatternInvDiagA.resize(numDofsFaceReduced, numDofsFaceReduced);

        // evaluate the acutal pattern
        for (int i = 0; i < invDiagA.N(); ++i)
        {
             occupationPatternInvDiagA.add(i, i);
        }

        occupationPatternInvDiagA.exportIdx(invDiagA);
    }

    /*!
     * \brief Resizes the  matrix ones and sets the matrix' sparsity pattern.
     */
    void setOnesPattern_(FaceToFaceMatrixBlock& ones, std::size_t numDofsFaceReduced)
    {
        // set the size of the sub-matrizes
        ones.setSize(numDofsFaceReduced, numDofsFaceReduced);

        // set occupation pattern of the coefficient matrix
        Dune::MatrixIndexSet occupationPattern;
        occupationPattern.resize(numDofsFaceReduced, numDofsFaceReduced);

        // evaluate the acutal pattern
        for (int i = 0; i < ones.N(); ++i)
        {
             occupationPattern.add(i, i);
        }

        occupationPattern.exportIdx(ones);
    }
};

} // end namespace Dumux

#endif
