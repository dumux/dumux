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
 * \brief Newton controller for multidomain problems
 */
#ifndef DUMUX_MULTIDOMAIN_NEWTON_CONTROLLER_HH
#define DUMUX_MULTIDOMAIN_NEWTON_CONTROLLER_HH

#include <dumux/nonlinear/newtoncontroller.hh>
#include "multidomainconvergencewriter.hh"

namespace Dumux
{
template <class TypeTag>
class MultiDomainNewtonController;

template <class TypeTag>
struct MultiDomainConvergenceWriter;

namespace Properties
{
// set default values for Newton for multidomain problems
// they can be overwritten in the parameter file
SET_INT_PROP(MultiDomain, NewtonTargetSteps, 8);
SET_INT_PROP(MultiDomain, NewtonMaxSteps, 15);
SET_SCALAR_PROP(MultiDomain, NewtonMaxRelativeShift, 1e-5);
}


/*!
 * \ingroup Newton
 * \ingroup MultidomainModel
 * \brief Reference implementation of a newton controller for coupled problems.
 *
 * If you want to specialize only some methods but are happy with
 * the defaults of the reference controller, derive your
 * controller from this class and simply overload the required
 * methods.
 */
template <class TypeTag>
class MultiDomainNewtonController : public NewtonController<TypeTag>
{
    typedef NewtonController<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonController) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, NewtonMethod) NewtonMethod;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, SplitAndMerge) SplitAndMerge;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, SubDomain1TypeTag) SubDomain1TypeTag;
    typedef typename GET_PROP_TYPE(TypeTag, SubDomain2TypeTag) SubDomain2TypeTag;

    typedef typename GET_PROP_TYPE(SubDomain1TypeTag, GridView) GridView1;
    typedef typename GET_PROP_TYPE(SubDomain2TypeTag, GridView) GridView2;

    typedef MultiDomainConvergenceWriter<TypeTag>  ConvergenceWriter;
    typedef typename GET_PROP_TYPE(TypeTag, LinearSolver) LinearSolver;

public:
    /*!
     * \brief Constructor
     *
     * \param problem The problem
     */
    MultiDomainNewtonController(const Problem &problem)
        : ParentType(problem)
		, endIterMsgStream_(std::ostringstream::out)
        , linearSolver_(problem)
        , convergenceWriter_(asImp_())
    {
//		  Writes out, where the relative tolerance is defined
        std::cout << "ParameterNewtonRelTol= "
        		<< PROP_DIAGNOSTIC(TypeTag, NewtonMaxRelativeShift)
        		<< ", "
        		<< GET_PROP_VALUE(TypeTag, NewtonMaxRelativeShift)
        		<< std::endl;
    }

    //! \copydoc ParentType::newtonUpdateShift()
    void newtonUpdateShift(const SolutionVector &uLastIter,
                           const SolutionVector &deltaU)
    {
        // calculate the relative error as the maximum relative
        // deflection in any degree of freedom.
        this->shift_ = 0;

        SolutionVector uNewI = uLastIter;
        uNewI -= deltaU;

        for (unsigned int i = 0; i < uLastIter.base().size(); ++i) {
            for (unsigned int j = 0; j < uLastIter.base()[i].size(); ++j) {
                Scalar vertexError = std::abs(deltaU.base()[i][j]);
                vertexError /= std::max<Scalar>(1.0, std::abs(uLastIter.base()[i][j] + uNewI.base()[i][j])/2);

                this->shift_ = std::max(this->shift_, vertexError);
            }
        }
    }

    void newtonUpdateRelError(const SolutionVector &uLastIter,
                              const SolutionVector &deltaU)
    DUNE_DEPRECATED_MSG("use newtonUpdateShift instead")
    { newtonUpdateShift(uLastIter, deltaU); }

    /*!
     * \brief Solve the linear system of equations
     *        \f$ \mathbf{A} x - b = 0\f$.
     *
     * \param A Coefficient matrix A
     * \param x Vector of unknowns
     * \param b Right hand side
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
            int converged = linearSolver_.solve(A.base(), x.base(), b.base());

#if HAVE_MPI
            // make sure all processes converged
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
#if HAVE_MPI
            // make sure all processes converged
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
            ms << e.what() << "M=" << A.base()[e.r][e.c];
            p.message(ms.str());
            throw p;
        }
        catch (const Dune::Exception &e) {
#if HAVE_MPI
            // make sure all processes converged
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
    * \param uCurrentIter The solution of the current iteration
    * \param uLastIter The solution of the last iteration
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

        newtonUpdateShift(uLastIter, deltaU);

        uCurrentIter = uLastIter;
        uCurrentIter -= deltaU;
    }

    /*!
    * \brief Indicates that one newton iteration was finished.
    *
    * \param uCurrentIter The solution of the current iteration
    * \param uLastIter The solution of the last iteration
    *
    */
    void newtonEndStep(SolutionVector &uCurrentIter, SolutionVector &uLastIter)
    {
        SplitAndMerge::splitSolVector(this->model_().curSol(),
                                      this->model_().sdModel1().curSol(),
                                      this->model_().sdModel2().curSol());

        ParentType::newtonEndStep(uCurrentIter, uLastIter);
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
            convergenceWriter_.beginIteration(sdGridView1_(), sdGridView2_());
            convergenceWriter_.writeFields(uLastIter, deltaU);
            convergenceWriter_.endIteration();
        }
    }

    /*!
     * \brief the subdomain gridviews
     */
    const GridView1 sdGridView1_() const
    { return this->problem_().sdGridView1(); }
    const GridView2 sdGridView2_() const
    { return this->problem_().sdGridView2(); }


private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    bool verbose_;

    std::ostringstream endIterMsgStream_;
    NewtonMethod *method_;

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

#endif // DUMUX_MULTIDOMAIN_NEWTON_CONTROLLER_HH
