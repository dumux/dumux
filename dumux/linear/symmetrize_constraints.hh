// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Linear
 * \brief Helper function to symmetrize row-constraints in constrained linear systems
 */
#ifndef DUMUX_LINEAR_SYMMETRIZE_CONSTRAINTS_HH
#define DUMUX_LINEAR_SYMMETRIZE_CONSTRAINTS_HH

#include <iostream>

#include <dune/common/hybridutilities.hh>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/multitypeblockmatrix.hh>

namespace Dumux::Detail {

// loop over dense matrix block: end of recursion (actually do some work)
template<class K, int rows, int cols, class VectorRow, class VectorCol>
void symmetrizeConstraintsImpl(Dune::FieldMatrix<K, rows, cols>& A,
                               VectorRow& bRow, const VectorCol& bCol,
                               const VectorCol& constrainedRows,
                               bool isDiagonal)
{
    for (int j = 0; j < A.M(); ++j)
    {
        // all rows in this column can be eliminated by a Dirichlet row
        const auto& constrainedRowsJ = constrainedRows[j];
        if (constrainedRowsJ > 0.5)
        {
            for (int i = 0; i < A.N(); ++i)
            {
                auto& Aij = A[i][j];
                auto& bi = bRow[i];
                const auto& bj = bCol[j];
                if (!((i == j) && isDiagonal))
                {
                    bi -= Aij*bj;
                    Aij = 0.0;
                }
            }
        }
    }
}

// recursively loop over sub-blocks
template<class MBlock, class VectorRow, class VectorCol>
void symmetrizeConstraintsImpl(Dune::BCRSMatrix<MBlock>& A,
                               VectorRow& bRow, const VectorCol& bCol,
                               const VectorCol& constrainedRows,
                               bool isDiagonal)
{
    const auto rowEnd = A.end();
    for (auto row = A.begin(); row != rowEnd; ++row)
    {
        const auto colEnd = row->end();
        for (auto col = row->begin(); col != colEnd; ++col)
        {
            const auto i = row.index();
            const auto j = col.index();

            auto& Aij = A[i][j];
            auto& bi = bRow[i];
            const auto& bj = bCol[j];
            const auto& dj = constrainedRows[j];

            symmetrizeConstraintsImpl(Aij, bi, bj, dj, ((i == j) && isDiagonal));
        }
    }
}

// recursively loop over sub-blocks
template<class... MBlock, class VectorRow, class VectorCol>
void symmetrizeConstraintsImpl(Dune::MultiTypeBlockMatrix<MBlock...>& A,
                               VectorRow& bRow, const VectorCol& bCol,
                               const VectorCol& constrainedRows,
                               bool isDiagonal)
{
    using namespace Dune::Hybrid;
    forEach(std::make_index_sequence<Dune::MultiTypeBlockMatrix<MBlock...>::N()>{}, [&](const auto i)
    {
        forEach(std::make_index_sequence<Dune::MultiTypeBlockMatrix<MBlock...>::M()>{}, [&](const auto j)
        {
            auto& Aij = A[i][j];
            auto& bi = bRow[i];
            const auto& bj = bCol[j];
            const auto& dj = constrainedRows[j];

            symmetrizeConstraintsImpl(Aij, bi, bj, dj, ((i == j) && isDiagonal));
        });
    });
}

} // end namespace Dumux::Detail

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief Symmetrize the constrained system Ax=b
 * \note The function addresses the case where some rows in the linear system have been replaced
 *       to fix a degree of freedom to a fixed value. A typical use case is strongly enforced Dirichlet constraints.
 *       To symmetrize the constraints in operator A, this function eliminates the corresponding column entries by
 *       bringing them to the right-hand side. Particularly, if matrix operator A was symmetric before incorporating
 *       constraints, calling this function restores the symmetry of operator A.
 * \note the constrainedRows vector is the same shape and type as b and contains 1 for each row
 *       that corresponds to constraints (e.g. strongly enforced Dirichlet dofs) and 0 elsewhere
 * \note this is a specialization for Dune::BCRSMatrix
 */
template<class MBlock, class Vector>
void symmetrizeConstraints(Dune::BCRSMatrix<MBlock>& A, Vector& b, const Vector& constrainedRows)
{
    Detail::symmetrizeConstraintsImpl(A, b, b, constrainedRows, true);
}

/*!
 * \ingroup Linear
 * \brief Symmetrize the constrained system Ax=b
 * \note The function addresses the case where some rows in the linear system have been replaced
 *       to fix a degree of freedom to a fixed value. A typical use case is strongly enforced Dirichlet constraints.
 *       To symmetrize the constraints in operator A, this function eliminates the corresponding column entries by
 *       bringing them to the right-hand side. Particularly, if matrix operator A was symmetric before incorporating
 *       constraints, calling this function restores the symmetry of operator A.
 * \note the constrainedRows vector is the same shape and type as b and contains 1 for each row
 *       that corresponds to constraints (e.g. strongly enforced Dirichlet dofs) and 0 elsewhere
 * \note this is a specialization for Dune::MultiTypeBlockMatrix
 */
template<class... MBlock, class Vector>
void symmetrizeConstraints(Dune::MultiTypeBlockMatrix<MBlock...>& A, Vector& b, const Vector& constrainedRows)
{
    Detail::symmetrizeConstraintsImpl(A, b, b, constrainedRows, true);
}

} // end namespace Dumux

#endif
