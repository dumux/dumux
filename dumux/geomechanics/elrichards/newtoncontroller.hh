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
 */
#ifndef DUMUX_ELRICHARDS_NEWTON_CONTROLLER_HH
#define DUMUX_ELRICHARDS_NEWTON_CONTROLLER_HH

#include <dumux/nonlinear/newtoncontroller.hh>

namespace Dumux {

/*!
* \brief An elrichards specific controller for the newton solver.
*
* This controller 'knows' what a 'physically meaningful' solution is
* which allows the newton method to abort quicker if the solution is
* way out of bounds.
 */
template <class TypeTag>
class ElRichardsNewtonController : public NewtonController<TypeTag>
{
    typedef NewtonController<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, NewtonController) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonMethod) NewtonMethod;

    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, LinearSolver) LinearSolver;

public:
    /*!
     * \brief Destructor
     */
    ElRichardsNewtonController(const Problem &problem)
    : ParentType(problem),linearSolver_(problem)
    {
        this->setTargetSteps(9);
        this->setMaxSteps(18);
    };

    void newtonUpdateRelError(const SolutionVector &uOld,
                              const SolutionVector &deltaU)
    {
        // calculate the relative error as the maximum relative
        // deflection in any degree of freedom.
        this->shift_ = 0;

        for (int i = 0; i < int(uOld.base().size()); ++i) {
            Scalar vertErr = std::abs(deltaU.base()[i]/(1.0 + std::abs((uOld.base()[i]) + uOld.base()[i] - deltaU.base()[i])/2));
            this->shift_ = std::max(this->shift_, vertErr);
        }

        this->shift_ = this->gridView_().comm().max(this->shift_);
    }

    void newtonUpdate(SolutionVector &uCurrentIter,
            const SolutionVector &uLastIter,
            const SolutionVector &deltaU)
    {
//        this->writeConvergence_(uLastIter, deltaU);

        newtonUpdateRelError(uLastIter, deltaU);

        uCurrentIter = uLastIter;
        uCurrentIter -= deltaU;

//        printvector(std::cout, deltaU, "new solution", "row", 12, 1, 3);
    }

    /*!
     * \brief Solve the linear system of equations \f$\mathbf{A}x - b = 0\f$.
     *
     * Throws NumericalProblem if the linear solver didn't
     * converge.
     *
     * \param A The matrix of the linear system of equations
     * \param x The vector which solves the linear system
     * \param b The right hand side of the linear system
     */
    void newtonSolveLinear(JacobianMatrix &A,
                           SolutionVector &x,
                           SolutionVector &b)
    {
        try {
            if (this->numSteps_ == 0)
            {
                Scalar norm2 = b.base().two_norm2();
                if (this->gridView_().comm().size() > 1)
                    norm2 = this->gridView_().comm().sum(norm2);

                initialAbsoluteError_ = std::sqrt(norm2);
                lastAbsoluteError_ = initialAbsoluteError_;
            }

            int converged = linearSolver_.solve(A.base(), x.base(), b.base());
//            printvector(std::cout, x.base(), "x", "row", 5, 1, 5);
//            printvector(std::cout, b.base(), "rhs", "row", 5, 1, 5);
//            Dune::writeMatrixToMatlab(A.base(), "matrix.txt");

            // make sure all processes converged
            int convergedRemote = converged;
            if (this->gridView_().comm().size() > 1)
                convergedRemote = this->gridView_().comm().min(converged);

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
            if (this->gridView_().comm().size() > 1)
                converged = this->gridView_().comm().min(converged);

            NumericalProblem p;
            std::string msg;
            std::ostringstream ms(msg);
            ms << e.what() << "M=" << A.base()[e.r][e.c];
            p.message(ms.str());
            throw p;
        }
        catch (const Dune::Exception &e) {
            // make sure all processes converged
            int converged = 0;
            if (this->gridView_().comm().size() > 1)
                converged = this->gridView_().comm().min(converged);

            NumericalProblem p;
            p.message(e.what());
            throw p;
        }
    }

    // absolute errors and tolerance
    Scalar absoluteError_;
    Scalar lastAbsoluteError_;
    Scalar initialAbsoluteError_;
    Scalar absoluteTolerance_;

    // the linear solver
    LinearSolver linearSolver_;

};
}

#endif
