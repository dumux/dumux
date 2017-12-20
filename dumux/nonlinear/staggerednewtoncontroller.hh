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
 * \ingroup Nonlinear
 * \ingroup StaggeredDiscretization
 * \brief A newton controller for staggered finite volume schemes
 */
#ifndef DUMUX_STAGGERED_NEWTON_CONTROLLER_HH
#define DUMUX_STAGGERED_NEWTON_CONTROLLER_HH

#include <dumux/common/properties.hh>
#include <dumux/common/exceptions.hh>

#include <dumux/nonlinear/newtoncontroller.hh>
#include <dumux/linear/linearsolveracceptsmultitypematrix.hh>
#include <dumux/linear/matrixconverter.hh>

namespace Dumux {

/*!
 * \ingroup Nonlinear
 * \ingroup StaggeredDiscretization
 * \brief A newton controller for staggered finite volume schemes
 */

template <class TypeTag>
class StaggeredNewtonController : public NewtonController<TypeTag>
{
    using ParentType = NewtonController<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    using GridView =  typename GET_PROP_TYPE(TypeTag, GridView);
    using Communicator = typename GridView::CollectiveCommunication;

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

    enum {
        numEqCellCenter = GET_PROP_VALUE(TypeTag, NumEqCellCenter),
        numEqFace = GET_PROP_VALUE(TypeTag, NumEqFace)
    };

public:
    using ParentType::ParentType;

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
    template<class LinearSolver, class JacobianMatrix, class SolutionVector>
    void solveLinearSystem(LinearSolver& ls,
                           JacobianMatrix& A,
                           SolutionVector& x,
                           SolutionVector& b)
    {
        try
        {
            if (this->numSteps_ == 0)
                this->initialResidual_ = b.two_norm();

            // check matrix sizes
            assert(A[cellCenterIdx][cellCenterIdx].N() == A[cellCenterIdx][faceIdx].N());
            assert(A[faceIdx][cellCenterIdx].N() == A[faceIdx][faceIdx].N());

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
            VectorConverter<SolutionVector>::retrieveValues(x, y);

            if (!converged)
                DUNE_THROW(NumericalProblem, "Linear solver did not converge");
        }
        catch (const Dune::Exception &e)
        {
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
    * \param assembler The assembler for Jacobian and residual
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
        if (this->enableShiftCriterion_)
            this->newtonUpdateShift(uLastIter, deltaU);

        if (this->useLineSearch_)
        {
            this->lineSearchUpdate_(assembler, uCurrentIter, uLastIter, deltaU);
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
                this->residualNorm_ = assembler.residualNorm(uCurrentIter);
                this->reduction_ = this->residualNorm_;
                this->reduction_ /= this->initialResidual_;
            }
        }

        // update the variables class to the new solution
        assembler.gridVariables().update(uCurrentIter);
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

            Scalar shiftAtDof = this->relativeShiftAtDof_(uLastIter[cellCenterIdx][i],
                                                            uNewI);
            this->shift_ = std::max(this->shift_, shiftAtDof);
        }
        for (int i = 0; i < int(uLastIter[faceIdx].size()); ++i) {
            auto uNewI = uLastIter[faceIdx][i];
            uNewI -= deltaU[faceIdx][i];

            Scalar shiftAtDof = this->relativeShiftAtDof_(uLastIter[faceIdx][i],
                                                            uNewI);
            this->shift_ = std::max(this->shift_, shiftAtDof);
        }

        if (this->communicator().size() > 1)
            this->shift_ = this->communicator().max(this->shift_);
    }

};

} // end namespace Dumux

#endif
