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
#include <dumux/linear/linearsolveracceptsmultitypematrix.hh>
#include "newtonconvergencewriter.hh"

#include "matrixconverter.hh"

#include <dune/common/tupleutility.hh>

namespace Dumux {

namespace Properties
{
    SET_PROP(StaggeredModel, LinearSolverBlockSize)
    {
        // LinearSolverAcceptsMultiTypeMatrix<T>::value
        // TODO: make somehow dependend? or only relevant for direct solvers?
    public:
        static constexpr auto value = 1;
    };
}
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

    enum {
        numEqCellCenter = GET_PROP_VALUE(TypeTag, NumEqCellCenter),
        numEqFace = GET_PROP_VALUE(TypeTag, NumEqFace)
    };

public:
    StaggeredNewtonController(const Problem &problem)
        : ParentType(problem)
    {}


    auto createMatrix(JacobianMatrix &A)
    {
        // copy the matrix and the vector to types the IterativeSolverBackend can handle

        using MatrixBlock = typename Dune::FieldMatrix<Scalar, 1, 1>;
        using SparseMatrix = typename Dune::BCRSMatrix<MatrixBlock>;

        // problem-agnostic version
        std::size_t numRows = 0;
        Dune::Hybrid::forEach(A[Dune::Indices::_0], [&numRows](const auto& subMatrixInFirstRow)
        {
            std::cout << "in matrix of first row:\n";
            std::cout << "numRows: " << subMatrixInFirstRow.N() << ", numCols: " << subMatrixInFirstRow.M() << std::endl;

            // The number of rows of the individual submatrice's block equals the respective number of equations.
            const auto numEq = std::decay_t<decltype(subMatrixInFirstRow)>::block_type::rows; // TODO: is this allways correct?
            numRows += numEq * subMatrixInFirstRow.M();
        });

        // /*const std::size_t*/ numRows = this->model_().numCellCenterDofs()*numEqCellCenter + this->model_().numFaceDofs()*numEqFace;


        // set the rowsizes
        Dune::MatrixIndexSet occupationPattern;
        occupationPattern.resize(numRows, numRows);

        // lambda function to fill the occupation pattern
        auto addIndices = [&occupationPattern](const auto& subMatrix, const std::size_t startRow, const std::size_t startCol)
        {
            using BlockType = typename std::decay_t<decltype(subMatrix)>::block_type;
            const auto blockSizeI = BlockType::rows;
            const auto blockSizeJ = BlockType::cols;
            for(auto row = subMatrix.begin(); row != subMatrix.end(); ++row)
                for(auto col = row->begin(); col != row->end(); ++col)
                    for(std::size_t i = 0; i < blockSizeI; ++i)
                        for(std::size_t j = 0; j < blockSizeJ; ++j)
                            occupationPattern.add(startRow + row.index()*blockSizeI + i, startCol + col.index()*blockSizeJ + j);

        };


        int rowIndex = 0;

        Dune::Hybrid::forEach(A, [&rowIndex, &addIndices, numRows](const auto& rowOfMultiTypeMatrix)
        {
            int colIndex = 0;

            std::cout << "processing one row of multitypeblockmatrix " << std::endl;
            Dune::Hybrid::forEach(rowOfMultiTypeMatrix, [&colIndex, &addIndices, &rowIndex, numRows](const auto& subMatrix)
            {
                std::cout << "numRows: " << subMatrix.N() << ", numCols: " << subMatrix.M() << std::endl;
                addIndices(subMatrix, rowIndex, colIndex);

                const auto numEq = std::decay_t<decltype(subMatrix)>::block_type::rows;
                colIndex += numEq*subMatrix.M();

                // if we arrived at the right side of the matrix, increase the row index
                if(colIndex == numRows)
                    rowIndex += numEq*subMatrix.N();
            });
        });

        auto M = SparseMatrix(numRows, numRows, SparseMatrix::random);
        occupationPattern.exportIdx(M);

        return M;
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
     * \param A The matrix of the linear system of equations
     * \param x The vector which solves the linear system
     * \param b The right hand side of the linear system
     */
    template<typename T = TypeTag>
    typename std::enable_if<!LinearSolverAcceptsMultiTypeMatrix<T>::value, void>::type
    newtonSolveLinear(JacobianMatrix &A,
                      SolutionVector &x,
                      SolutionVector &b)
    {
        // createMatrix(A);
        try
        {
            if (this->numSteps_ == 0)
                this->initialResidual_ = b.two_norm();

            // copy the matrix and the vector to types the IterativeSolverBackend can handle
            using MatrixBlock = typename Dune::FieldMatrix<Scalar, 1, 1>;
            using SparseMatrix = typename Dune::BCRSMatrix<MatrixBlock>;

            // get the new matrix sizes
            std::size_t numRows = numEqCellCenter*A[cellCenterIdx][cellCenterIdx].N() + numEqFace*A[faceIdx][cellCenterIdx].N();
            std::size_t numCols = numEqCellCenter*A[cellCenterIdx][cellCenterIdx].M() + numEqFace*A[cellCenterIdx][faceIdx].M();

            // check matrix sizes
            assert(A[cellCenterIdx][cellCenterIdx].N() == A[cellCenterIdx][faceIdx].N());
            assert(A[faceIdx][cellCenterIdx].N() == A[faceIdx][faceIdx].N());
            assert(numRows == numCols);

            // create the bcrs matrix the IterativeSolver backend can handle
            const auto M = MatrixConverter<JacobianMatrix>::multiTypeToBCRSMatrix(A);

            // create the vector the IterativeSolver backend can handle
            using VectorBlock = typename Dune::FieldVector<Scalar, 1>;
            using BlockVector = typename Dune::BlockVector<VectorBlock>;

            BlockVector y, bTmp;
            y.resize(numRows);
            bTmp.resize(numCols);
            for (std::size_t i = 0; i < b[cellCenterIdx].N(); ++i)
                for (std::size_t j = 0; j < numEqCellCenter; ++j)
                    bTmp[i*numEqCellCenter + j] = b[cellCenterIdx][i][j];
            for (std::size_t i = 0; i < b[faceIdx].N(); ++i)
                for (std::size_t j = 0; j < numEqFace; ++j)
                    bTmp[i*numEqFace + j + b[cellCenterIdx].N()*numEqCellCenter] = b[faceIdx][i][j];

            // printmatrix(std::cout, M, "", "");

            // solve
            const bool converged = this->linearSolver_.solve(M, y, bTmp);

            // copy back the result y into x
            for (std::size_t i = 0; i < x[cellCenterIdx].N(); ++i)
                for (std::size_t j = 0; j < numEqCellCenter; ++j)
                    x[cellCenterIdx][i][j] = y[i*numEqCellCenter + j];
            for (std::size_t i = 0; i < x[faceIdx].N(); ++i)
                for (std::size_t j = 0; j < numEqFace; ++j)
                    x[faceIdx][i][j] = y[i*numEqFace + j + x[cellCenterIdx].N()*numEqCellCenter];

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
     * \brief Solve the linear system of equations \f$\mathbf{A}x - b = 0\f$.
     *
     * Throws Dumux::NumericalProblem if the linear solver didn't
     * converge.
     *
     * \param A The matrix of the linear system of equations
     * \param x The vector which solves the linear system
     * \param b The right hand side of the linear system
     */
    template<typename T = TypeTag>
    typename std::enable_if<LinearSolverAcceptsMultiTypeMatrix<T>::value, void>::type
    newtonSolveLinear(JacobianMatrix &A,
                      SolutionVector &x,
                      SolutionVector &b)
    {
        try
        {
            if (this->numSteps_ == 0)
                this->initialResidual_ = b.two_norm();

            bool converged = this->linearSolver_.solve(A, x, b);

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
