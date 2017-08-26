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
 * \brief Reference implementation of a controller class for the Newton solver.
 *
 * Usually this controller should be sufficient.
 */
#ifndef DUMUX_NEWTON_CONTROLLER_HH
#define DUMUX_NEWTON_CONTROLLER_HH

#include <dumux/common/basicproperties.hh>
#include <dumux/common/propertysystem.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/math.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/linear/seqsolverbackend.hh>

#include "newtonmethod.hh"

namespace Dumux
{
template <class TypeTag>
class NewtonController;

namespace Properties
{
//! Specifies the implementation of the Newton controller
NEW_PROP_TAG(NewtonController);

//! Specifies whether the Jacobian matrix should only be reassembled
//! if the current solution deviates too much from the evaluation point
NEW_PROP_TAG(ImplicitEnablePartialReassemble);

//! Specify whether the jacobian matrix of the last iteration of a
//! time step should be re-used as the jacobian of the first iteration
//! of the next time step.
NEW_PROP_TAG(ImplicitEnableJacobianRecycling);

/*!
 * \brief Specifies whether the update should be done using the line search
 *        method instead of the plain Newton method.
 *
 * Whether this property has any effect depends on whether the line
 * search method is implemented for the actual model's Newton
 * controller's update() method. By default line search is not used.
 */
NEW_PROP_TAG(NewtonUseLineSearch);

//! indicate whether the absolute residual should be used in the residual criterion
NEW_PROP_TAG(NewtonEnableAbsoluteResidualCriterion);

//! indicate whether the shift criterion should be used
NEW_PROP_TAG(NewtonEnableShiftCriterion);

//! the value for the maximum relative shift below which convergence is declared
NEW_PROP_TAG(NewtonMaxRelativeShift);

//! the value for the maximum absolute residual below which convergence is declared
NEW_PROP_TAG(NewtonMaxAbsoluteResidual);

//! indicate whether the residual criterion should be used
NEW_PROP_TAG(NewtonEnableResidualCriterion);

//! the value for the residual reduction below which convergence is declared
NEW_PROP_TAG(NewtonResidualReduction);

//! indicate whether both of the criteria should be satisfied to declare convergence
NEW_PROP_TAG(NewtonSatisfyResidualAndShiftCriterion);

//! indicate whether after each newton step the solution is chopped, e.g. to
//! physically sensible values (if a chopper is implemented for this model)
NEW_PROP_TAG(NewtonEnableChop);

/*!
 * \brief The number of iterations at which the Newton method
 *        should aim at.
 *
 * This is used to control the time-step size. The heuristic used
 * is to scale the last time-step size by the deviation of the
 * number of iterations used from the target steps.
 */
NEW_PROP_TAG(NewtonTargetSteps);

//! Number of maximum iterations for the Newton method.
NEW_PROP_TAG(NewtonMaxSteps);

//! set default values
SET_TYPE_PROP(NewtonMethod, NewtonController, NewtonController<TypeTag>);
SET_BOOL_PROP(NewtonMethod, NewtonUseLineSearch, false);
SET_BOOL_PROP(NewtonMethod, NewtonEnableAbsoluteResidualCriterion, false);
SET_BOOL_PROP(NewtonMethod, NewtonEnableShiftCriterion, true);
SET_BOOL_PROP(NewtonMethod, NewtonEnableResidualCriterion, false);
SET_BOOL_PROP(NewtonMethod, NewtonSatisfyResidualAndShiftCriterion, false);
SET_BOOL_PROP(NewtonMethod, NewtonEnableChop, false);
SET_SCALAR_PROP(NewtonMethod, NewtonMaxRelativeShift, 1e-8);
SET_SCALAR_PROP(NewtonMethod, NewtonMaxAbsoluteResidual, 1e-5);
SET_SCALAR_PROP(NewtonMethod, NewtonResidualReduction, 1e-5);
SET_INT_PROP(NewtonMethod, NewtonTargetSteps, 10);
SET_INT_PROP(NewtonMethod, NewtonMaxSteps, 18);

} // end namespace Properties

/*!
 * \ingroup Newton
 * \brief A reference implementation of a Newton controller specific
 *        for the box scheme.
 *
 * If you want to specialize only some methods but are happy with the
 * defaults of the reference controller, derive your controller from
 * this class and simply overload the required methods.
 */
template <class TypeTag>
class NewtonController
{
    using Scalar =  typename GET_PROP_TYPE(TypeTag, Scalar);
    using Implementation =  typename GET_PROP_TYPE(TypeTag, NewtonController);
    using GridView =  typename GET_PROP_TYPE(TypeTag, GridView);
    using Communicator = typename GridView::CollectiveCommunication;
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);

    static constexpr int numEq = GET_PROP_VALUE(TypeTag, NumEq);

public:
    /*!
     * \brief Constructor for stationary problems
     */
    NewtonController(const Communicator& comm)
    : comm_(comm)
    , endIterMsgStream_(std::ostringstream::out)
    {
        initParams_();
    }

    /*!
     * \brief Constructor for stationary problems
     */
    NewtonController(const Communicator& comm, std::shared_ptr<TimeLoop<Scalar>> timeLoop)
    : comm_(comm)
    , timeLoop_(timeLoop)
    , endIterMsgStream_(std::ostringstream::out)
    {
        initParams_();
    }

    const Communicator& communicator() const
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
     * \brief Returns true if another iteration should be done.
     *
     * \param uCurrentIter The solution of the current Newton iteration
     */
    template<class SolutionVector>
    bool newtonProceed(const SolutionVector &uCurrentIter)
    {
        if (numSteps_ < 2)
            return true; // we always do at least two iterations
        else if (asImp_().newtonConverged()) {
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
     * \brief Returns true if the error of the solution is below the
     *        tolerance.
     */
    bool newtonConverged() const
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
     * \brief Called before the Newton method is applied to an
     *        non-linear system of equations.
     *
     * \param u The initial solution
     */
    template<class SolutionVector>
    void newtonBegin(const SolutionVector &u)
    {
        numSteps_ = 0;
    }

    /*!
     * \brief Indicates the beginning of a Newton iteration.
     */
    void newtonBeginStep()
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
     * \brief Returns the number of steps done since newtonBegin() was
     *        called.
     */
    int newtonNumSteps()
    { return numSteps_; }

    /*!
     * \brief Update the maximum relative shift of the solution compared to
     *        the previous iteration.
     *
     * \param uLastIter The current iterative solution
     * \param deltaU The difference between the current and the next solution
     */
    template<class SolutionVector>
    void newtonUpdateShift(const SolutionVector &uLastIter,
                           const SolutionVector &deltaU)
    {
        shift_ = 0;

        for (int i = 0; i < int(uLastIter.size()); ++i) {
            typename SolutionVector::block_type uNewI = uLastIter[i];
            uNewI -= deltaU[i];

            Scalar shiftAtDof = relativeShiftAtDof_(uLastIter[i], uNewI);
            using std::max;
            shift_ = max(shift_, shiftAtDof);
        }

        if (communicator().size() > 1)
            shift_ = communicator().max(shift_);
    }

    /*!
     * \brief Assemble the linear system of equations \f$\mathbf{A}x - b = 0\f$.
     *
     * \param assembler The jacobian assembler
     */
    template<class JacobianAssembler, class SolutionVector>
    void assembleLinearSystem(JacobianAssembler& assembler,
                              const SolutionVector& uCurrentIter)
    {
        assembler.assembleJacobianAndResidual(uCurrentIter);
    }

    /*!
     * \brief Solve the linear system of equations \f$\mathbf{A}x - b = 0\f$.
     *
     * Throws NumericalProblem if the linear solver didn't
     * converge.
     *
     * \param ls The linear solver to be used
     * \param A The matrix of the linear system of equations
     * \param x The vector which solves the linear system
     * \param b The right hand side of the linear system
     */
    template<class LinearSolver, class JacobianMatrix, class SolutionVector>
    void solveLinearSystem(LinearSolver& ls,
                           JacobianMatrix& A,
                           SolutionVector& x,
                           SolutionVector& b)
    {
        try {
            if (numSteps_ == 0)
            {
                Scalar norm2 = b.two_norm2();
                if (communicator().size() > 1)
                    norm2 = communicator().sum(norm2);

                using std::sqrt;
                initialResidual_ = sqrt(norm2);
            }

            //! Copy into a standard block vector. This is necessary for all model _not_ using a FieldVector<Scalar, numEq> as
            //! primary variables vector in combination with UMFPack or SuperLU as their interfaces are hard coded
            //! to this field vector type in Dune ISTL
            //! Could be avoided for vectors that already have the right type using SFINAE
            //! but it shouldn't impact performance too much
            Dune::BlockVector<NumEqVector> xTmp; xTmp.resize(b.size());
            Dune::BlockVector<NumEqVector> bTmp(xTmp);
            for (int i = 0; i < b.size(); ++i)
                for (int j = 0; j < numEq; ++j)
                    bTmp[i][j] = b[i][j];

            int converged = ls.solve(A, xTmp, bTmp);

            for (int i = 0; i < x.size(); ++i)
                for (int j = 0; j < numEq; ++j)
                    x[i][j] = xTmp[i][j];

            // make sure all processes converged
            int convergedRemote = converged;
            if (communicator().size() > 1)
                convergedRemote = communicator().min(converged);

            if (!converged) {
                DUNE_THROW(NumericalProblem,
                           "Linear solver did not converge");
            }
            else if (!convergedRemote) {
                DUNE_THROW(NumericalProblem,
                           "Linear solver did not converge on a remote process");
            }
        }
        catch (Dune::MatrixBlockError e) {
            // make sure all processes converged
            int converged = 0;
            if (communicator().size() > 1)
                converged = communicator().min(converged);

            NumericalProblem p;
            std::string msg;
            std::ostringstream ms(msg);
            ms << e.what() << "M=" << A[e.r][e.c];
            p.message(ms.str());
            throw p;
        }
        catch (const Dune::Exception &e) {
            // make sure all processes converged
            int converged = 0;
            if (communicator().size() > 1)
                converged = communicator().min(converged);

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
    template<class JacobianAssembler, class SolutionVector>
    void newtonUpdate(JacobianAssembler& assembler,
                      SolutionVector &uCurrentIter,
                      const SolutionVector &uLastIter,
                      const SolutionVector &deltaU)
    {
        if (enableShiftCriterion_)
            newtonUpdateShift(uLastIter, deltaU);

        if (useLineSearch_)
        {
            lineSearchUpdate_(assembler, uCurrentIter, uLastIter, deltaU);
        }
        else {
            for (unsigned int i = 0; i < uLastIter.size(); ++i) {
                uCurrentIter[i] = uLastIter[i];
                uCurrentIter[i] -= deltaU[i];
            }

            if (enableResidualCriterion_)
            {
                residualNorm_ = assembler.globalResidual(uCurrentIter);
                reduction_ = residualNorm_;
                reduction_ /= initialResidual_;
            }
        }

        // update the variables class to the new solution
        assembler.gridVariables().update(uCurrentIter);
    }

    /*!
     * \brief Indicates that one Newton iteration was finished.
     *
     * \param assembler The jacobian assembler
     * \param uCurrentIter The solution after the current Newton iteration
     * \param uLastIter The solution at the beginning of the current Newton iteration
     */
    template<class JacobianAssembler, class SolutionVector>
    void newtonEndStep(JacobianAssembler& assembler,
                       SolutionVector &uCurrentIter,
                       const SolutionVector &uLastIter)
    {
        ++numSteps_;

        if (verbose())
        {
            std::cout << "\rNewton iteration " << numSteps_ << " done";
            if (enableShiftCriterion_)
                std::cout << ", maximum relative shift = " << shift_;
            if (enableResidualCriterion_ && enableAbsoluteResidualCriterion_)
                std::cout << ", residual = " << residualNorm_;
            else if (enableResidualCriterion_)
                std::cout << ", residual reduction = " << reduction_;
            std::cout << endIterMsg().str() << "\n";
        }
        endIterMsgStream_.str("");

        // When the Newton iterations are done: ask the model to check whether it makes sense
        // TODO: how do we realize this? -> do this here in the newton controller
        // model_().checkPlausibility();
    }

    /*!
     * \brief Called if the Newton method broke down.
     *
     * This method is called _after_ newtonEnd()
     */
    template<class Assembler, class SolutionVector>
    void newtonFail(Assembler& assembler, SolutionVector& u)
    {
        if (!assembler.localResidual().isStationary())
        {
            // set solution to previous solution
            u = assembler.prevSol();

            // reset the grid variables to the previous solution
            assembler.gridVariables().resetTimeStep(u);

            std::cout << "Newton solver did not converge with dt = "
                      << timeLoop_->timeStepSize() << " seconds. Retrying with time step of "
                      << timeLoop_->timeStepSize()/2 << " seconds\n";

            // try again with dt = dt/2
            timeLoop_->setTimeStepSize(timeLoop_->timeStepSize()/2);
        }
        else
            DUNE_THROW(Dune::MathError, "Newton solver did not converge");
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

    std::ostringstream &endIterMsg()
    { return endIterMsgStream_; }

    /*!
     * \brief Specifies if the Newton method ought to be chatty.
     */
    void setVerbose(bool val)
    { verbose_ = val; }

    /*!
     * \brief Returns true if the Newton method ought to be chatty.
     */
    bool verbose() const
    { return verbose_ && communicator().rank() == 0; }

protected:

    // returns the actual implementation for the controller we do
    // it this way in order to allow "poor man's virtual methods",
    // i.e. methods of subclasses which can be called by the base
    // class.
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    void initParams_()
    {
        enablePartialReassemble_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Implicit, EnablePartialReassemble);
        enableJacobianRecycling_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Implicit, EnableJacobianRecycling);

        useLineSearch_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Newton, UseLineSearch);
        enableAbsoluteResidualCriterion_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Newton, EnableAbsoluteResidualCriterion);
        enableShiftCriterion_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Newton, EnableShiftCriterion);
        enableResidualCriterion_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Newton, EnableResidualCriterion)
                                   || enableAbsoluteResidualCriterion_;
        satisfyResidualAndShiftCriterion_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Newton, SatisfyResidualAndShiftCriterion);
        if (!enableShiftCriterion_ && !enableResidualCriterion_)
        {
            DUNE_THROW(Dune::NotImplemented,
                       "at least one of NewtonEnableShiftCriterion or "
                       << "NewtonEnableResidualCriterion has to be set to true");
        }

        setMaxRelativeShift(GET_PARAM_FROM_GROUP(TypeTag, Scalar, Newton, MaxRelativeShift));
        setMaxAbsoluteResidual(GET_PARAM_FROM_GROUP(TypeTag, Scalar, Newton, MaxAbsoluteResidual));
        setResidualReduction(GET_PARAM_FROM_GROUP(TypeTag, Scalar, Newton, ResidualReduction));
        setTargetSteps(GET_PARAM_FROM_GROUP(TypeTag, int, Newton, TargetSteps));
        setMaxSteps(GET_PARAM_FROM_GROUP(TypeTag, int, Newton, MaxSteps));

        verbose_ = true;
        numSteps_ = 0;
    }

    template<class JacobianAssembler, class SolutionVector>
    void lineSearchUpdate_(const JacobianAssembler& assembler,
                           SolutionVector &uCurrentIter,
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

            residualNorm_ = assembler.globalResidual(uCurrentIter);
            reduction_ = residualNorm_;
            reduction_ /= initialResidual_;

            if (reduction_ < lastReduction_ || lambda <= 0.125) {
                endIterMsg() << ", residual reduction " << lastReduction_ << "->"  << reduction_ << "@lambda=" << lambda;
                return;
            }

            // try with a smaller update
            lambda /= 2.0;
        }
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
                               const PrimaryVariables &priVars2)
    {
        Scalar result = 0.0;
        using std::abs;
        using std::max;
        for (int j = 0; j < numEq; ++j) {
            Scalar eqErr = abs(priVars1[j] - priVars2[j]);
            eqErr /= max<Scalar>(1.0,abs(priVars1[j] + priVars2[j])/2);

            result = max(result, eqErr);
        }
        return result;
    }

    // The grid view's communicator
    const Communicator& comm_;

    std::shared_ptr<TimeLoop<Scalar>> timeLoop_;

    // message stream to be displayed at the end of iterations
    std::ostringstream endIterMsgStream_;

    // switches on/off verbosity
    bool verbose_;

    // shift criterion variables
    Scalar shift_;
    Scalar lastShift_;
    Scalar shiftTolerance_;

    // residual criterion variables
    Scalar reduction_;
    Scalar residualNorm_;
    Scalar lastReduction_;
    Scalar initialResidual_;
    Scalar reductionTolerance_;
    Scalar residualTolerance_;

    // optimal number of iterations we want to achieve
    int targetSteps_;
    // maximum number of iterations we do before giving up
    int maxSteps_;
    // actual number of steps done so far
    int numSteps_;

    // further parameters
    bool enablePartialReassemble_;
    bool enableJacobianRecycling_;
    bool useLineSearch_;
    bool enableAbsoluteResidualCriterion_;
    bool enableShiftCriterion_;
    bool enableResidualCriterion_;
    bool satisfyResidualAndShiftCriterion_;
};

} // namespace Dumux

#endif
