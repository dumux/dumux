// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
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
* \brief Reference implementation of a newton controller for coupled problems.
*/
#ifndef DUMUX_MULTIDOMAIN_NEWTON_CONTROLLER_HH
#define DUMUX_MULTIDOMAIN_NEWTON_CONTROLLER_HH

#include <dumux/common/exceptions.hh>
#include <dumux/linear/linearsolverproperties.hh>
#include <dumux/linear/boxlinearsolver.hh>

#include <dune/istl/overlappingschwarz.hh>
#include <dune/istl/schwarz.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include "dune/istl/owneroverlapcopy.hh"

#include <dune/istl/io.hh>
#include <dune/common/mpihelper.hh>
#include <iostream>
#include <boost/format.hpp>

#include <dune/pdelab/backend/istlsolverbackend.hh>

#include "multidomainconvergencewriter.hh"

/*!
 * \file
 * \brief Additional properties required for the coupled Newton controller
 */
namespace Dumux
{
template <class TypeTag>
class MultiDomainNewtonController;

namespace Properties
{

//! Specifies the implementation of the Newton controller
NEW_PROP_TAG(NewtonController);

//! Specifies the type of the actual Newton method
NEW_PROP_TAG(NewtonMethod);

//! specifies whether the convergence rate and the global residual
//! gets written out to disk for every newton iteration (default is false)
NEW_PROP_TAG(NewtonWriteConvergence);

/*!
 * \brief Specifies whether the update should be done using the line search
 *        method instead of the plain Newton method.
 *
 * Whether this property has any effect depends on wether the line
 * search method is implemented for the actual model's Newton
 * controller's update() method. By default line search is not used.
 */
NEW_PROP_TAG(NewtonUseLineSearch);

//! the value for the relative error below which convergence is declared
NEW_PROP_TAG(NewtonRelTolerance);

/*!
 * \brief The number of iterations at which the Newton method
 *        should aim at.
 *
 * This is used to control the time step size. The heuristic used
 * is to scale the last time step size by the deviation of the
 * number of iterations used from the target steps.
 */
NEW_PROP_TAG(NewtonTargetSteps);

//! Number of maximum iterations for the Newton method.
NEW_PROP_TAG(NewtonMaxSteps);

// set default values for Newton
// they can be overwritten in the parameter file
SET_INT_PROP(MultiDomain, NewtonTargetSteps, 8);
SET_INT_PROP(MultiDomain, NewtonMaxSteps, 15);
SET_SCALAR_PROP(MultiDomain, NewtonRelTolerance, 1e-5);
SET_BOOL_PROP(MultiDomain, NewtonWriteConvergence, false);
}


/*!
 * \brief Reference implementation of a newton controller for coupled problems.
 *
 * If you want to specialize only some methods but are happy with
 * the defaults of the reference controller, derive your
 * controller from this class and simply overload the required
 * methods.
 */
template <class TypeTag>
class MultiDomainNewtonController
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonController) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonMethod) NewtonMethod;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, SubDomain1TypeTag) SubTypeTag1;
    typedef typename GET_PROP_TYPE(TypeTag, SubDomain2TypeTag) SubTypeTag2;

    typedef typename GET_PROP_TYPE(SubTypeTag1, GridView) GridView1;
    typedef typename GET_PROP_TYPE(SubTypeTag2, GridView) GridView2;

    typedef MultiDomainConvergenceWriter<TypeTag>  ConvergenceWriter;
    typedef typename GET_PROP_TYPE(TypeTag, LinearSolver) LinearSolver;

public:
    MultiDomainNewtonController(const Problem &problem)
        : endIterMsgStream_(std::ostringstream::out)
        , linearSolver_(problem)
        , convergenceWriter_(asImp_())
    {
        verbose_ = true;
        numSteps_ = 0;

        // maximum tolerated deflection between two iterations
        setRelTolerance(GET_PARAM_FROM_GROUP(TypeTag, Scalar, Newton, RelTolerance));
        setTargetSteps(GET_PARAM_FROM_GROUP(TypeTag, Scalar, Newton, TargetSteps));
        setMaxSteps(GET_PARAM_FROM_GROUP(TypeTag, Scalar, Newton, MaxSteps));

//		  Writes out, where the relative tolerance is defined
//        std::cout << "ParameterNewtonRelTol= " << PROP_DIAGNOSTIC(TypeTag, NewtonRelTolerance) << std::endl;
    };

    /*!
     * \brief Destructor
     */
    ~MultiDomainNewtonController()
    {
    };

    /*!
     * \brief Specifies if the newton method ought to be chatty.
     *
     * \param val docme
     *
     */
    void setVerbose(bool val)
    { verbose_ = val; }

    /*!
     * \brief Returns true if the newton method ought to be chatty.
     */
    bool verbose() const
    { return verbose_ && gridView_().comm().rank() == 0; }

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
//    void setAbsTolerance(Scalar tolerance)
//    { absoluteTolerance_ = tolerance; }

    /*!
     * \brief Set the number of iterations at which the Newton method
     *        should aim at.
     *
     * This is used to control the time step size. The heuristic used
     * is to scale the last time step size by the deviation of the
     * number of iterations used from the target steps.
     *
     * \param targetSteps Number of iterations which are considered "optimal"
     */
    void setTargetSteps(int targetSteps)
    { targetSteps_ = targetSteps; }

    /*!
     * \brief Set the number of iterations after which the Newton
     *        method gives up.
     */
    void setMaxSteps(int maxSteps)
    { maxSteps_ = maxSteps; }

    /*!
    * \brief Returns true if another iteration should be done.
     *
     * \param uCurrentIter The solution of the current newton iteration
    */
    bool newtonProceed(SolutionVector &uCurrentIter)
    {
        if (numSteps_ < 2)
            return true; // we always do at least two iterations
        else
            if (numSteps_ >= maxSteps_)
            {
                // we have exceeded the allowed number of steps.  if the
                // relative error was reduced by a factor of at least 5,
                // we proceed anyway
                return error_*4.0 < lastError_ && !asImp_().newtonConverged();
            }
            else
                if (asImp_().newtonConverged())
                    return false; // we are below the desired tolerance

        return true; // do another round
    }

    /*!
    * \brief Returns true if the error of the solution is below the
    *        tolerance.
    */
    bool newtonConverged() const
    { return error_ <= tolerance_; }

    /*!
    * \brief Called before the newton method is applied to an
    *        non-linear system of equations.
    *
    * \param method docme
    * \param uCurrentIter docme
    *
    */
    void newtonBegin(NewtonMethod &method, const SolutionVector &uCurrentIter)
    {
        method_ = &method;
        numSteps_ = 0;

        if (GET_PARAM_FROM_GROUP(TypeTag, bool, Newton, WriteConvergence))
            convergenceWriter_.beginTimestep();
    }

    /*!
    * \brief Indidicates the beginning of a newton iteration.
    */
    void newtonBeginStep()
    { lastError_ = error_; }

    /*!
    * \brief Returns the number of steps done since newtonBegin() was
    *        called.
    */
    int newtonNumSteps()
    { return numSteps_; }


    /*!
    * \brief Update the error of the solution compared to the
    *        previous iteration.
    *
    * \param uLastIter docme
    * \param deltaU docme
    *
    */
    void newtonUpdateRelError(const SolutionVector &uLastIter,
                              const SolutionVector &deltaU)
    {
        // calculate the relative error as the maximum relative
        // deflection in any degree of freedom.
        error_ = 0;

        SolutionVector uNewI = uLastIter;
        uNewI -= deltaU;

        for (unsigned int i = 0; i < uLastIter.size(); ++i) {
            for (unsigned int j = 0; j < uLastIter[i].size(); ++j) {
                Scalar vertexError = std::abs(deltaU[i][j]);
                vertexError /= std::max<Scalar>(1.0, std::abs(uLastIter[i][j] + uNewI[i][j])/2);

                error_ = std::max(error_, vertexError);
            }
        }
    }


    /*!
     * \brief Solve the linear system of equations \f$ \mathbf{A}x - b
     *        = 0\f$.
     *
     * \param A docme
     * \param x docme
     * \param b docme
     *
     * Throws Dumux::NumericalProblem if the linear solver didn't
     * converge.
     */
    template <class Matrix, class Vector>
    void newtonSolveLinear(Matrix &A,
                           Vector &x,
                           Vector &b)
    {
        // if the deflection of the newton method is large, we do not
        // need to solve the linear approximation accurately. Assuming
        // that the initial value for the delta vector u is quite
        // close to the final value, a reduction of 6 orders of
        // magnitude in the defect should be sufficient...
        try {
            int converged = linearSolver_.solve(A, x, b);

            // make sure all processes converged
#if HAVE_MPI
            int convergedSend = 1;
            MPI_Allreduce(/*sendBuf=*/&convergedSend,
                          /*recvBuf=*/&converged,
                          /*count=*/1,
                          MPI_INT,
                          MPI_MIN,
                          MPI_COMM_WORLD);
#endif
            if (!converged) {
                DUNE_THROW(NumericalProblem,
                           "Linear solver did not converge");
            }
        }
        catch (const Dune::MatrixBlockError &e) {
            // make sure all processes converged
#if HAVE_MPI
            int convergedSend = 0;
            int converged;

            MPI_Allreduce(/*sendBuf=*/&convergedSend,
                          /*recvBuf=*/&converged,
                          /*count=*/1,
                          MPI_INT,
                          MPI_MIN,
                          MPI_COMM_WORLD);
#endif

            Dumux::NumericalProblem p;
            std::string msg;
            std::ostringstream ms(msg);
            ms << e.what() << "M=" << A[e.r][e.c];
            p.message(ms.str());
            throw p;
        }
        catch (const Dune::Exception &e) {
            // make sure all processes converged
#if HAVE_MPI
            int convergedSend = 0;
            int converged;

            MPI_Allreduce(/*sendBuf=*/&convergedSend,
                          /*recvBuf=*/&converged,
                          /*count=*/1,
                          MPI_INT,
                          MPI_MIN,
                          MPI_COMM_WORLD);
#endif

            Dumux::NumericalProblem p;
            p.message(e.what());
            throw p;
        }
    }

    /*!
    * \brief Update the current solution function with a delta vector.
    *
    * The error estimates required for the newtonConverged() and
    * newtonProceed() methods should be updated here.
    *
    * Different update strategies, such as line search and chopped
    * updates can be implemented. The default behaviour is just to
    * subtract deltaU from uLastIter.
    *
    * \param uCurrentIter docme
    * \param uLastIter   The solution of the last iteration
    * \param deltaU The delta as calculated from solving the linear
    *               system of equations. This parameter also stores
    *               the updated solution.
    *
    */
    void newtonUpdate(SolutionVector &uCurrentIter,
                      const SolutionVector &uLastIter,
                      const SolutionVector &deltaU)
    {
        if (GET_PARAM_FROM_GROUP(TypeTag, bool, Newton, WriteConvergence)) {
            writeConvergence_(uLastIter, deltaU);
        }

        newtonUpdateRelError(uLastIter, deltaU);

        uCurrentIter = uLastIter;
        uCurrentIter -= deltaU;
    }

    /*!
    * \brief Indicates that one newton iteration was finished.
    *
    * \param uCurrentIter docme
    * \param uLastIter The solution of the last iteration
    *
    */
    void newtonEndStep(SolutionVector &uCurrentIter, SolutionVector &uLastIter)
    {
        typedef Dumux::SplitAndMerge<TypeTag> Common;

        Common::splitSolVector(this->model().curSol(),
                               this->model().subModel1().curSol(),
                               this->model().subModel2().curSol());

        ++numSteps_;

        if (verbose())
            std::cout << "\rNewton iteration " << numSteps_ << " done: "
                      << "error=" << error_ << endIterMsg().str() << "\n";
        endIterMsgStream_.str("");
    }

    /*!
    * \brief Indicates that we're done solving the non-linear system of equations.
    */
    void newtonEnd()
    {
        if (GET_PARAM_FROM_GROUP(TypeTag, bool, Newton, WriteConvergence))
            convergenceWriter_.endTimestep();
    }

    /*!
    * \brief Called if the newton method broke down.
    *
    * This method is called _after_ newtonEnd()
    */
    void newtonFail()
    {
        numSteps_ = std::max(maxSteps_, targetSteps_*2);
    }

    /*!
    * \brief Called when the newton method was sucessful.
    *
    * This method is called _after_ newtonEnd()
    */
    void newtonSucceed()
    {
    }

    /*!
    * \brief Suggest a new time stepsize based on the old time step size.
    *
    * The default behaviour is to suggest the old time step size
    * scaled by the ratio between the target iterations and the
    * iterations required to actually solve the last time step.
    *
    * \param oldTimeStep docme
    *
    */
    Scalar suggestTimeStepSize(Scalar oldTimeStep) const
    {
        // be agressive reducing the timestep size but
        // conservative when increasing it. the rationale is
        // that we want to avoid failing in the next newton
        // iteration which would require another linearization
        // of the problem.
        if (numSteps_ > targetSteps_) {
            Scalar percent = ((Scalar) numSteps_ - targetSteps_)/targetSteps_;
            return oldTimeStep/(1.0 + percent);
        }
        else {
            /*Scalar percent = (Scalar(1))/targetSteps_;
              return oldTimeStep*(1 + percent);
             */
            Scalar percent = ((Scalar) targetSteps_ - numSteps_)/targetSteps_;
            return oldTimeStep*(1.0 + percent/1.2);
        }
    }

    /*!
    * \brief Returns a reference to the current newton method
    *        which is controlled by this controller.
    */
    NewtonMethod &method()
    { return *method_; }

    /*!
    * \brief Returns a reference to the current newton method
    *        which is controlled by this controller.
    */
    const NewtonMethod &method() const
    { return *method_; }

    /*!
    * \brief Returns a reference to the current numeric model.
    */
    Model &model()
    { return method_->model(); }

    /*!
    * \brief Returns a const reference to the current numeric model.
    */
    const Model &model() const
    { return method_->model(); }

    std::ostringstream &endIterMsg()
    { return endIterMsgStream_; }

    const GridView1 &gridView1() const
    { return problem_().gridView1(); }
    const GridView2 &gridView2() const
    { return problem_().gridView2(); }

    /*!
    * \brief the convergence writer produces the output
    *
    * \param uLastIter The solution of the last iteration
    * \param deltaU The delta as calculated from solving the linear
    *               system of equations. This parameter also stores
    *               the updated solution.
    *
    */
    void writeConvergence_(const SolutionVector &uLastIter,
                           const SolutionVector &deltaU)
    {
        if (GET_PARAM_FROM_GROUP(TypeTag, bool, Newton, WriteConvergence)) {
            convergenceWriter_.beginIteration(gridView1_(), gridView2_());
            convergenceWriter_.writeFields(uLastIter, deltaU);
            convergenceWriter_.endIteration();
        };
    };

    /*!
     * \brief Returns a copy of the the grid view.
     */
    const GridView gridView_() const
    { return problem_().gridView(); }

    /*!
     * \brief the subdomain gridviews
     */
    const GridView1 gridView1_() const
    { return problem_().subProblem1().gridView(); }
    const GridView2 gridView2_() const
    { return problem_().subProblem2().gridView(); }

    /*!
     * \brief the coupled problem
     */
    Problem &problem_()
        { return method_->problem(); }
    const Problem &problem_() const
        { return method_->problem(); }

    // returns the actual implementation for the controller we do
    // it this way in order to allow "poor man's virtual methods",
    // i.e. methods of subclasses which can be called by the base
    // class.
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    bool verbose_;

    std::ostringstream endIterMsgStream_;
    NewtonMethod *method_;

    Scalar tolerance_;
    Scalar error_;
    Scalar lastError_;

    // optimal number of iterations we want to achieve
    int targetSteps_;
    // maximum number of iterations we do before giving up
    int maxSteps_;
    // actual number of steps done so far
    int numSteps_;

    // the linear solver
    LinearSolver linearSolver_;

    ConvergenceWriter convergenceWriter_;
};

} // namespace Dumux

#endif
