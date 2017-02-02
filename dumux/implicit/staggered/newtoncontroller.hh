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
 * \brief A 2p1cni specific controller for the newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
#ifndef DUMUX_STAGGERED_NEWTON_CONTROLLER_HH
#define DUMUX_STAGGERED_NEWTON_CONTROLLER_HH

#include "properties.hh"

#include <dumux/nonlinear/newtoncontroller.hh>
#include "newtonconvergencewriter.hh"

namespace Dumux {
/*!
 * \ingroup PNMModel
 * \brief A PNM specific controller for the newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
template <class TypeTag>
class StaggeredNewtonController : public NewtonController<TypeTag>
{
    typedef NewtonController<TypeTag> ParentType;
    typedef NewtonConvergenceWriter<TypeTag> StaggeredNewtonConvergenceWriter;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

public:
    StaggeredNewtonController(const Problem &problem)
        : ParentType(problem)
    {}

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
    void newtonSolveLinear(JacobianMatrix &A,
                           SolutionVector &x,
                           SolutionVector &b)
    {
        try {
            if (this->numSteps_ == 0)
            {
                Scalar norm2 = b.two_norm2();
                if (this->gridView_().comm().size() > 1)
                    norm2 = this->gridView_().comm().sum(norm2);

                this->initialResidual_ = std::sqrt(norm2);
            }

            int converged = this->linearSolver_.solve(A, x, b);

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

            Dumux::NumericalProblem p;
            std::string msg;
            std::ostringstream ms(msg);
//             ms << e.what() << "M=" << A[e.r][e.c];
//             p.message(ms.str());
//             throw p;
        }
        catch (const Dune::Exception &e) {
            // make sure all processes converged
            int converged = 0;
            if (this->gridView_().comm().size() > 1)
                converged = this->gridView_().comm().min(converged);

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
        if (this->enableShiftCriterion_)
            this->newtonUpdateShift(uLastIter, deltaU);

        this->writeConvergence_(uLastIter, deltaU);

        if (this->useLineSearch_)
        {
            this->lineSearchUpdate_(uCurrentIter, uLastIter, deltaU);
        }
        else {
            for (unsigned int i = 0; i < uLastIter[cellCenterIdx].size(); ++i) {
                uCurrentIter[cellCenterIdx][i] = uLastIter[cellCenterIdx][i];
                uCurrentIter[cellCenterIdx][i] -= deltaU[cellCenterIdx][i];
            }
            for (unsigned int i = 0; i < uLastIter[faceIdx].size(); ++i) {
                uCurrentIter[faceIdx][i] = uLastIter[faceIdx][i];
                uCurrentIter[faceIdx][i] -= deltaU[faceIdx][i];
            }

            if (this->enableResidualCriterion_)
            {
                SolutionVector tmp(uLastIter);
                this->reduction_ = this->method().model().globalResidual(tmp, uCurrentIter);
                this->reduction_ /= this->initialResidual_;
            }
        }
    }

     /*!
     * \brief Update the maximum relative shift of the solution compared to
     *        the previous iteration.
     *
     * \param uLastIter The current iterative solution
     * \param deltaU The difference between the current and the next solution
     */
    void newtonUpdateShift(const SolutionVector &uLastIter,
                           const SolutionVector &deltaU)
    {
        this->shift_ = 0;

        for (int i = 0; i < int(uLastIter[cellCenterIdx].size()); ++i) {
            auto uNewI = uLastIter[cellCenterIdx][i];
            uNewI -= deltaU[cellCenterIdx][i];

            Scalar shiftAtDof = this->model_().relativeShiftAtDof(uLastIter[cellCenterIdx][i],
                                                            uNewI);
            this->shift_ = std::max(this->shift_, shiftAtDof);
        }
        for (int i = 0; i < int(uLastIter[faceIdx].size()); ++i) {
            auto uNewI = uLastIter[faceIdx][i];
            uNewI -= deltaU[faceIdx][i];

            Scalar shiftAtDof = this->model_().relativeShiftAtDof(uLastIter[faceIdx][i],
                                                            uNewI);
            this->shift_ = std::max(this->shift_, shiftAtDof);
        }

        if (this->gridView_().comm().size() > 1)
            this->shift_ = this->gridView_().comm().max(this->shift_);
    }

};
}

#endif
