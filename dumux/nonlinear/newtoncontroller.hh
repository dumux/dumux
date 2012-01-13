// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/****************************************************************************
 *   Copyright (C) 2008-2011 by Andreas Lauser                               *
 *   Copyright (C) 2008-2010 by Bernd Flemisch                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
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

#include <dumux/common/propertysystem.hh>
#include <dumux/io/vtkmultiwriter.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/math.hh>
#include <dumux/io/vtkmultiwriter.hh>
#include <dumux/linear/boxlinearsolver.hh>

#include "newtonconvergencewriter.hh"

namespace Dumux
{
template <class TypeTag>
class NewtonController;

template <class TypeTag>
class NewtonConvergenceWriter;

namespace Properties
{
//! Specifies the implementation of the Newton controller
NEW_PROP_TAG(NewtonController);

//! Specifies the type of the actual Newton method
NEW_PROP_TAG(NewtonMethod);

//! Specifies the type of a solution
NEW_PROP_TAG(SolutionVector);

//! Specifies the type of a vector of primary variables at a degree of freedom
NEW_PROP_TAG(PrimaryVariables);

//! Specifies the type of a global Jacobian matrix
NEW_PROP_TAG(JacobianMatrix);

//! Specifies the type of the Vertex mapper
NEW_PROP_TAG(VertexMapper);

//! specifies the type of the time manager
NEW_PROP_TAG(TimeManager);

//! specifies whether the convergence rate and the global residual
//! gets written out to disk for every Newton iteration (default is false)
NEW_PROP_TAG(NewtonWriteConvergence);

//! Specifies whether the Jacobian matrix should only be reassembled
//! if the current solution deviates too much from the evaluation point
NEW_PROP_TAG(EnablePartialReassemble);

/*!
 * \brief Specifies whether the update should be done using the line search
 *        method instead of the plain Newton method.
 *
 * Whether this property has any effect depends on whether the line
 * search method is implemented for the actual model's Newton
 * controller's update() method. By default line search is not used.
 */
NEW_PROP_TAG(NewtonUseLineSearch);

//! indicate whether the relative error should be used
NEW_PROP_TAG(NewtonEnableRelativeCriterion);

//! the value for the relative error below which convergence is declared
NEW_PROP_TAG(NewtonRelTolerance);

//! indicate whether the absolute error should be used
NEW_PROP_TAG(NewtonEnableAbsoluteCriterion);

//! the value for the absolute error reduction below which convergence is declared
NEW_PROP_TAG(NewtonAbsTolerance);

//! indicate whether both of the criteria should be satisfied to declare convergence
NEW_PROP_TAG(NewtonSatisfyAbsAndRel);

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

// set default values
SET_TYPE_PROP(NewtonMethod, NewtonController, Dumux::NewtonController<TypeTag>);
SET_BOOL_PROP(NewtonMethod, NewtonWriteConvergence, false);
SET_BOOL_PROP(NewtonMethod, NewtonUseLineSearch, false);
SET_BOOL_PROP(NewtonMethod, NewtonEnableRelativeCriterion, true);
SET_BOOL_PROP(NewtonMethod, NewtonEnableAbsoluteCriterion, false);
SET_BOOL_PROP(NewtonMethod, NewtonSatisfyAbsAndRel, false);
SET_SCALAR_PROP(NewtonMethod, NewtonRelTolerance, 1e-8);
SET_SCALAR_PROP(NewtonMethod, NewtonAbsTolerance, 1e-5);
SET_INT_PROP(NewtonMethod, NewtonTargetSteps, 10);
SET_INT_PROP(NewtonMethod, NewtonMaxSteps, 18);

}

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
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonController) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonMethod) NewtonMethod;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, VertexMapper) VertexMapper;

    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    typedef NewtonConvergenceWriter<TypeTag> ConvergenceWriter;

    typedef typename GET_PROP_TYPE(TypeTag, LinearSolver) LinearSolver;

public:
    /*!
     * \brief
     */
    NewtonController(const Problem &problem)
        : endIterMsgStream_(std::ostringstream::out)
        , convergenceWriter_(asImp_())
        , linearSolver_(problem)
    {
        enablePartialReassemble_ = GET_PARAM(TypeTag, bool, EnablePartialReassemble);
        enableJacobianRecycling_ = GET_PARAM(TypeTag, bool, EnableJacobianRecycling);

        useLineSearch_ = GET_PARAM(TypeTag, bool, Newton, UseLineSearch);
        enableRelativeCriterion_ = GET_PARAM(TypeTag, bool, Newton, EnableRelativeCriterion);
        enableAbsoluteCriterion_ = GET_PARAM(TypeTag, bool, Newton, EnableAbsoluteCriterion);
        satisfyAbsAndRel_ = GET_PARAM(TypeTag, bool, Newton, SatisfyAbsAndRel);
        if (!enableRelativeCriterion_ && !enableAbsoluteCriterion_)
        {
            DUNE_THROW(Dune::NotImplemented, "at least one of NewtonEnableRelativeCriterion or "
                    << "NewtonEnableAbsoluteCriterion has to be set to true");
        }

        setRelTolerance(GET_PARAM(TypeTag, Scalar, Newton, RelTolerance));
        setAbsTolerance(GET_PARAM(TypeTag, Scalar, Newton, AbsTolerance));
        setTargetSteps(GET_PARAM(TypeTag, int, Newton, TargetSteps));
        setMaxSteps(GET_PARAM(TypeTag, int, Newton, MaxSteps));

        verbose_ = true;
        numSteps_ = 0;
    };

    /*!
     * \brief Destructor
     */
    ~NewtonController()
    {
    };

    /*!
     * \brief Set the maximum acceptable difference for convergence of
     *        any primary variable between two iterations.
     *
     * \param tolerance The maximum relative error between two Newton
     *                  iterations at which the scheme is considered
     *                  finished
     */
    void setRelTolerance(Scalar tolerance)
    { tolerance_ = tolerance; }

    /*!
     * \brief Set the maximum acceptable residual norm reduction.
     *
     * \param tolerance The maximum reduction of the residual norm
     *                  at which the scheme is considered finished
     */
    void setAbsTolerance(Scalar tolerance)
    { absoluteTolerance_ = tolerance; }

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
    bool newtonProceed(const SolutionVector &uCurrentIter)
    {
        if (numSteps_ < 2)
            return true; // we always do at least two iterations
        else if (asImp_().newtonConverged()) {
            return false; // we are below the desired tolerance
        }
        else if (numSteps_ >= maxSteps_) {
            // we have exceeded the allowed number of steps.  if the
            // relative error was reduced by a factor of at least 4,
            // we proceed even if we are above the maximum number of
            // steps
            if (enableRelativeCriterion_)
                return error_*4.0 < lastError_;
            else
                return absoluteError_*4.0 < lastAbsoluteError_;
        }

        return true;
    }

    /*!
     * \brief Returns true if the error of the solution is below the
     *        tolerance.
     */
    bool newtonConverged() const
    {
        if (enableRelativeCriterion_ && !enableAbsoluteCriterion_)
        {
            return error_ <= tolerance_;
        }
        else if (!enableRelativeCriterion_ && enableAbsoluteCriterion_)
        {
            return absoluteError_ <= absoluteTolerance_;
        }
        else if (satisfyAbsAndRel_)
        {
            return error_ <= tolerance_
                    && absoluteError_ <= absoluteTolerance_;
        }
        else
        {
            return error_ <= tolerance_
                    || absoluteError_ <= absoluteTolerance_;
        }

        return false;
    }

    /*!
     * \brief Called before the Newton method is applied to an
     *        non-linear system of equations.
     *
     * \param method The object where the NewtonMethod is executed
     * \param u The initial solution
     */
    void newtonBegin(NewtonMethod &method, const SolutionVector &u)
    {
        method_ = &method;
        numSteps_ = 0;

        if (GET_PARAM(TypeTag, bool, Newton, WriteConvergence))
            convergenceWriter_.beginTimestep();
    }

    /*!
     * \brief Indicates the beginning of a Newton iteration.
     */
    void newtonBeginStep()
    {
        lastError_ = error_;
        lastAbsoluteError_ = absoluteError_;
    }

    /*!
     * \brief Returns the number of steps done since newtonBegin() was
     *        called.
     */
    int newtonNumSteps()
    { return numSteps_; }

    /*!
     * \brief Update the relative error of the solution compared to
     *        the previous iteration.
     *
     * The relative error can be seen as a norm of the difference
     * between the current and the next iteration.
     *
     * \param uLastIter The current iterative solution
     * \param deltaU The difference between the current and the next solution
     */
    void newtonUpdateRelError(const SolutionVector &uLastIter,
                              const SolutionVector &deltaU)
    {
        // calculate the relative error as the maximum relative
        // deflection in any degree of freedom.
        error_ = 0;

        for (int i = 0; i < int(uLastIter.size()); ++i) {
            PrimaryVariables uNewI = uLastIter[i];
            uNewI -= deltaU[i];

            Scalar vertError = model_().relativeErrorVertex(i,
                                                            uLastIter[i],
                                                            uNewI);
            error_ = std::max(error_, vertError);

        }

        if (gridView_().comm().size() > 1)
            error_ = gridView_().comm().max(error_);
    }

    /*!
     * \brief Solve the linear system of equations \f$\mathbf{A}x - b = 0\f$.
     *
     * Throws Dumux::NumericalProblem if the linear solver didn't
     * converge.
     *
     * \param A The matrix of the linear system of equations
     * \param x The vector which solves the linear system
     * \param b The right hand side of the linear system
     */
    void newtonSolveLinear(const JacobianMatrix &A,
                           SolutionVector &x,
                           const SolutionVector &b)
    {
        try {
            if (numSteps_ == 0)
            {
                Scalar norm2 = b.two_norm2();
                if (gridView_().comm().size() > 1)
                    norm2 = gridView_().comm().sum(norm2);

                initialAbsoluteError_ = std::sqrt(norm2);
                lastAbsoluteError_ = initialAbsoluteError_;
            }

            int converged = linearSolver_.solve(A, x, b);

            // make sure all processes converged
            int convergedRemote = converged;
            if (gridView_().comm().size() > 1)
                convergedRemote = gridView_().comm().min(converged);

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
            if (gridView_().comm().size() > 1)
                converged = gridView_().comm().min(converged);

            Dumux::NumericalProblem p;
            std::string msg;
            std::ostringstream ms(msg);
            ms << e.what() << "M=" << A[e.r][e.c];
            p.message(ms.str());
            throw p;
        }
        catch (const Dune::Exception &e) {
            // make sure all processes converged
            int converged = 0;
            if (gridView_().comm().size() > 1)
                converged = gridView_().comm().min(converged);

            Dumux::NumericalProblem p;
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
        writeConvergence_(uLastIter, deltaU);

        if (enableRelativeCriterion_ || enablePartialReassemble_)
            newtonUpdateRelError(uLastIter, deltaU);

        // compute the vertex and element colors for partial reassembly
        if (enablePartialReassemble_) {
            Scalar minReasmTol, tmp, reassembleTol;
            minReasmTol = 0.1*tolerance_;
            tmp = Dumux::geometricMean(error_, minReasmTol);
            reassembleTol = Dumux::geometricMean(error_, tmp);
            reassembleTol = std::max(reassembleTol, minReasmTol);
            this->model_().jacobianAssembler().updateDiscrepancy(uLastIter, deltaU);
            this->model_().jacobianAssembler().computeColors(reassembleTol);
        }

        if (useLineSearch_)
        {
            lineSearchUpdate_(uCurrentIter, uLastIter, deltaU);
        }
        else {
            uCurrentIter = uLastIter;
            uCurrentIter -= deltaU;

            if (enableAbsoluteCriterion_)
            {
                SolutionVector tmp(uLastIter);
                absoluteError_ = this->method().model().globalResidual(tmp, uCurrentIter);
                absoluteError_ /= initialAbsoluteError_;
            }
        }

    }

    /*!
     * \brief Indicates that one Newton iteration was finished.
     *
     * \param uCurrentIter The solution after the current Newton iteration
     * \param uLastIter The solution at the beginning of the current Newton iteration
     */
    void newtonEndStep(const SolutionVector &uCurrentIter,
                       const SolutionVector &uLastIter)
    {
        ++numSteps_;

        if (verbose())
        {
            std::cout << "\rNewton iteration " << numSteps_ << " done";
            if (enableRelativeCriterion_)
                std::cout << ", relative error = " << error_;
            if (enableAbsoluteCriterion_)
                std::cout << ", absolute error = " << absoluteError_;
            std::cout << endIterMsg().str() << "\n";
        }

        endIterMsgStream_.str("");
    }

    /*!
     * \brief Indicates that we're done solving the non-linear system
     *        of equations.
     */
    void newtonEnd()
    {
        if (GET_PARAM(TypeTag, bool, Newton, WriteConvergence))
            convergenceWriter_.endTimestep();
    }

    /*!
     * \brief Called if the Newton method broke down.
     *
     * This method is called _after_ newtonEnd()
     */
    void newtonFail()
    {
        model_().jacobianAssembler().reassembleAll();
        numSteps_ = targetSteps_*2;
    }

    /*!
     * \brief Called when the Newton method was successful.
     *
     * This method is called _after_ newtonEnd()
     */
    void newtonSucceed()
    {
        if (enableJacobianRecycling_)
            model_().jacobianAssembler().setMatrixReuseable(true);
        else
            model_().jacobianAssembler().reassembleAll();
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
        else {
            Scalar percent = Scalar(targetSteps_ - numSteps_)/targetSteps_;
            return oldTimeStep*(1.0 + percent/1.2);
        }
    }

    /*!
     * \brief Returns a reference to the current Newton method
     *        which is controlled by this controller.
     */
    NewtonMethod &method()
    { return *method_; }

    /*!
     * \brief Returns a reference to the current Newton method
     *        which is controlled by this controller.
     */
    const NewtonMethod &method() const
    { return *method_; }

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
    { return verbose_ && gridView_().comm().rank() == 0; }

protected:
    /*!
     * \brief Returns a reference to the grid view.
     */
    const GridView &gridView_() const
    { return problem_().gridView(); }

    /*!
     * \brief Returns a reference to the vertex mapper.
     */
    const VertexMapper &vertexMapper_() const
    { return model_().vertexMapper(); }

    /*!
     * \brief Returns a reference to the problem.
     */
    Problem &problem_()
    { return method_->problem(); }

    /*!
     * \brief Returns a reference to the problem.
     */
    const Problem &problem_() const
    { return method_->problem(); }

    /*!
     * \brief Returns a reference to the time manager.
     */
    TimeManager &timeManager_()
    { return problem_().timeManager(); }

    /*!
     * \brief Returns a reference to the time manager.
     */
    const TimeManager &timeManager_() const
    { return problem_().timeManager(); }

    /*!
     * \brief Returns a reference to the problem.
     */
    Model &model_()
    { return problem_().model(); }

    /*!
     * \brief Returns a reference to the problem.
     */
    const Model &model_() const
    { return problem_().model(); }

    // returns the actual implementation for the controller we do
    // it this way in order to allow "poor man's virtual methods",
    // i.e. methods of subclasses which can be called by the base
    // class.
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    void writeConvergence_(const SolutionVector &uLastIter,
                           const SolutionVector &deltaU)
    {
        if (GET_PARAM(TypeTag, bool, Newton, WriteConvergence)) {
            convergenceWriter_.beginIteration(this->gridView_());
            convergenceWriter_.writeFields(uLastIter, deltaU);
            convergenceWriter_.endIteration();
        }
    };

    void lineSearchUpdate_(SolutionVector &uCurrentIter,
                           const SolutionVector &uLastIter,
                           const SolutionVector &deltaU)
    {
       Scalar lambda = 1.0;
       SolutionVector tmp(uLastIter);

       while (true) {
           uCurrentIter = deltaU;
           uCurrentIter *= -lambda;
           uCurrentIter += uLastIter;

           // calculate the residual of the current solution
           absoluteError_ = this->method().model().globalResidual(tmp, uCurrentIter);
           absoluteError_ /= initialAbsoluteError_;

           if (absoluteError_ < lastAbsoluteError_ || lambda <= 0.125) {
               this->endIterMsg() << ", defect " << lastAbsoluteError_ << "->"  << absoluteError_ << "@lambda=" << lambda;
               return;
           }

           // try with a smaller update
           lambda /= 2.0;
       }
    }

    std::ostringstream endIterMsgStream_;

    NewtonMethod *method_;
    bool verbose_;

    ConvergenceWriter convergenceWriter_;

    // relative errors and tolerance
    Scalar error_;
    Scalar lastError_;
    Scalar tolerance_;

    // absolute errors and tolerance
    Scalar absoluteError_;
    Scalar lastAbsoluteError_;
    Scalar initialAbsoluteError_;
    Scalar absoluteTolerance_;

    // optimal number of iterations we want to achieve
    int targetSteps_;
    // maximum number of iterations we do before giving up
    int maxSteps_;
    // actual number of steps done so far
    int numSteps_;

    // the linear solver
    LinearSolver linearSolver_;

private:
    bool enablePartialReassemble_;
    bool enableJacobianRecycling_;
    bool useLineSearch_;
    bool enableRelativeCriterion_;
    bool enableAbsoluteCriterion_;
    bool satisfyAbsAndRel_;
};
} // namespace Dumux

#endif
