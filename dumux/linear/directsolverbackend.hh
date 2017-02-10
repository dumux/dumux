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
 * \brief Dumux multidimension direct solver backend
 */
#ifndef DUMUX_STAGGEREDGRID_DIRECT_SOLVER_BACKEND_HH
#define DUMUX_STAGGEREDGRID_DIRECT_SOLVER_BACKEND_HH

//#include <dune/istl/superlu.hh>
#include <dune/istl/umfpack.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/basicproperties.hh>
#include <dumux/linear/linearsolverproperties.hh>

namespace Dumux {

#if HAVE_UMFPACK
/*! \brief UMFPackBackend for MultiTypeMatrices
 *         Copies the coupled matrix into a single BCRSMatrix.
 *         Very slow!! Only wise to use for verification purposes.
 */
template <class TypeTag>
class StaggeredGridUMFPackBackend
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    enum {
        numEqCellCenter = GET_PROP_VALUE(TypeTag, NumEqCellCenter),
        numEqFace = GET_PROP_VALUE(TypeTag, NumEqFace)
    };

    typename GET_PROP(TypeTag, DofTypeIndices)::CellCenterIdx cellCenterIdx;
    typename GET_PROP(TypeTag, DofTypeIndices)::FaceIdx faceIdx;

public:

  StaggeredGridUMFPackBackend(const Problem& problem) {}

  template<class Matrix, class Vector>
  bool solve(const Matrix& A, Vector& x, const Vector& b)
  {
    // copy the matrix and the vector to types the UMFPackBackend can handle
    using MatrixBlock = typename Dune::FieldMatrix<Scalar, 1, 1>;
    using SparseMatrix = typename Dune::BCRSMatrix<MatrixBlock>;

    // get the new matrix sizes
    std::size_t numRows = numEqCellCenter*A[cellCenterIdx][cellCenterIdx].N() + numEqFace*A[faceIdx][cellCenterIdx].N();
    std::size_t numCols = numEqCellCenter*A[cellCenterIdx][cellCenterIdx].M() + numEqFace*A[cellCenterIdx][faceIdx].M();

    // check matrix sizes
    assert(A[cellCenterIdx][cellCenterIdx].N() == A[cellCenterIdx][faceIdx].N());
    assert(A[faceIdx][cellCenterIdx].N() == A[faceIdx][faceIdx].N());
    assert(numRows == numCols);

    // create the bcrs matrix the UMFPack backend can handle
    auto M = SparseMatrix(numRows, numCols, SparseMatrix::random);

    // set the rowsizes
    // A11 and A12
    for (auto row = A[cellCenterIdx][cellCenterIdx].begin(); row != A[cellCenterIdx][cellCenterIdx].end(); ++row)
        for (std::size_t i = 0; i < numEqCellCenter; ++i)
            M.setrowsize(numEqCellCenter*row.index() + i, row->size()*numEqCellCenter);
    for (auto row = A[cellCenterIdx][faceIdx].begin(); row != A[cellCenterIdx][faceIdx].end(); ++row)
        for (std::size_t i = 0; i < numEqCellCenter; ++i)
            M.setrowsize(numEqCellCenter*row.index() + i, M.getrowsize(numEqCellCenter*row.index() + i) + row->size()*numEqFace);
    // A21 and A22
    for (auto row = A[faceIdx][cellCenterIdx].begin(); row != A[faceIdx][cellCenterIdx].end(); ++row)
        for (std::size_t i = 0; i < numEqFace; ++i)
            M.setrowsize(numEqFace*row.index() + i + A[cellCenterIdx][cellCenterIdx].N()*numEqCellCenter, row->size()*numEqCellCenter);
    for (auto row = A[faceIdx][faceIdx].begin(); row != A[faceIdx][faceIdx].end(); ++row)
        for (std::size_t i = 0; i < numEqFace; ++i)
            M.setrowsize(numEqFace*row.index() + i + A[cellCenterIdx][cellCenterIdx].N()*numEqCellCenter, M.getrowsize(numEqFace*row.index() + i + A[cellCenterIdx][cellCenterIdx].N()*numEqCellCenter) + row->size()*numEqFace);
    M.endrowsizes();

    // set the indices
    for (auto row = A[cellCenterIdx][cellCenterIdx].begin(); row != A[cellCenterIdx][cellCenterIdx].end(); ++row)
        for (auto col = row->begin(); col != row->end(); ++col)
            for (std::size_t i = 0; i < numEqCellCenter; ++i)
                for (std::size_t j = 0; j < numEqCellCenter; ++j)
                    M.addindex(row.index()*numEqCellCenter + i, col.index()*numEqCellCenter + j);

    for (auto row = A[cellCenterIdx][faceIdx].begin(); row != A[cellCenterIdx][faceIdx].end(); ++row)
        for (auto col = row->begin(); col != row->end(); ++col)
            for (std::size_t i = 0; i < numEqCellCenter; ++i)
                for (std::size_t j = 0; j < numEqFace; ++j)
                    M.addindex(row.index()*numEqCellCenter + i, col.index()*numEqFace + j + A[cellCenterIdx][cellCenterIdx].M()*numEqCellCenter);

    for (auto row = A[faceIdx][cellCenterIdx].begin(); row != A[faceIdx][cellCenterIdx].end(); ++row)
        for (auto col = row->begin(); col != row->end(); ++col)
            for (std::size_t i = 0; i < numEqFace; ++i)
                for (std::size_t j = 0; j < numEqCellCenter; ++j)
                    M.addindex(row.index()*numEqFace + i + A[cellCenterIdx][cellCenterIdx].N()*numEqCellCenter, col.index()*numEqCellCenter + j);

    for (auto row = A[faceIdx][faceIdx].begin(); row != A[faceIdx][faceIdx].end(); ++row)
        for (auto col = row->begin(); col != row->end(); ++col)
            for (std::size_t i = 0; i < numEqFace; ++i)
                for (std::size_t j = 0; j < numEqFace; ++j)
                    M.addindex(row.index()*numEqFace + i + A[cellCenterIdx][cellCenterIdx].N()*numEqCellCenter, col.index()*numEqFace + j + A[cellCenterIdx][cellCenterIdx].M()*numEqCellCenter);
    M.endindices();

    // copy values
    for (auto row = A[cellCenterIdx][cellCenterIdx].begin(); row != A[cellCenterIdx][cellCenterIdx].end(); ++row)
        for (auto col = row->begin(); col != row->end(); ++col)
            for (std::size_t i = 0; i < numEqCellCenter; ++i)
                for (std::size_t j = 0; j < numEqCellCenter; ++j)
                    M[row.index()*numEqCellCenter + i][col.index()*numEqCellCenter + j] = A[cellCenterIdx][cellCenterIdx][row.index()][col.index()][i][j];

    for (auto row = A[cellCenterIdx][faceIdx].begin(); row != A[cellCenterIdx][faceIdx].end(); ++row)
        for (auto col = row->begin(); col != row->end(); ++col)
            for (std::size_t i = 0; i < numEqCellCenter; ++i)
                for (std::size_t j = 0; j < numEqFace; ++j)
                    M[row.index()*numEqCellCenter + i][col.index()*numEqFace + j + A[cellCenterIdx][cellCenterIdx].M()*numEqCellCenter] = A[cellCenterIdx][faceIdx][row.index()][col.index()][i][j];

    for (auto row = A[faceIdx][cellCenterIdx].begin(); row != A[faceIdx][cellCenterIdx].end(); ++row)
        for (auto col = row->begin(); col != row->end(); ++col)
            for (std::size_t i = 0; i < numEqFace; ++i)
                for (std::size_t j = 0; j < numEqCellCenter; ++j)
                    M[row.index()*numEqFace + i + A[cellCenterIdx][cellCenterIdx].N()*numEqCellCenter][col.index()*numEqCellCenter + j] = A[faceIdx][cellCenterIdx][row.index()][col.index()][i][j];

    for (auto row = A[faceIdx][faceIdx].begin(); row != A[faceIdx][faceIdx].end(); ++row)
        for (auto col = row->begin(); col != row->end(); ++col)
            for (std::size_t i = 0; i < numEqFace; ++i)
                for (std::size_t j = 0; j < numEqFace; ++j)
                    M[row.index()*numEqFace + i + A[cellCenterIdx][cellCenterIdx].N()*numEqCellCenter][col.index()*numEqFace + j + A[cellCenterIdx][cellCenterIdx].M()*numEqCellCenter] = A[faceIdx][faceIdx][row.index()][col.index()][i][j];

    // create the vector the UMFPack backend can handle
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

    int verbosity = GET_PARAM_FROM_GROUP(TypeTag, int, LinearSolver, Verbosity);
    Dune::UMFPack<SparseMatrix> solver(M, verbosity > 0);

    solver.apply(y, bTmp, result_);

    std::size_t size = y.size();
    for (std::size_t i = 0; i < size; i++)
    {
        if (std::isnan(y[i][0]) || std::isinf(y[i][0]))
        {
            result_.converged = false;
            break;
        }
    }

    // copy back the result y into x
    for (std::size_t i = 0; i < x[cellCenterIdx].N(); ++i)
        for (std::size_t j = 0; j < numEqCellCenter; ++j)
            x[cellCenterIdx][i][j] = y[i*numEqCellCenter + j];
    for (std::size_t i = 0; i < x[faceIdx].N(); ++i)
        for (std::size_t j = 0; j < numEqFace; ++j)
            x[faceIdx][i][j] = y[i*numEqFace + j + x[cellCenterIdx].N()*numEqCellCenter];

//     printmatrix(std::cout,M, "umfpack", "");

    return result_.converged;
  }

  const Dune::InverseOperatorResult& result() const
  {
    return result_;
  }

private:
  Dune::InverseOperatorResult result_;
};

#endif // HAVE_UMFPACK

} // end namespace Dumux

#endif
