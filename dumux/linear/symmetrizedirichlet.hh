// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup Linear
 * \brief Helper function to symmetrize strong Dirichlet conditions
 */
#ifndef DUMUX_LINEAR_SYMMETRIZE_DIRICHLET_HH
#define DUMUX_LINEAR_SYMMETRIZE_DIRICHLET_HH

#include <iostream>

#include <dune/common/hybridutilities.hh>
#include <dune/istl/blocklevel.hh>

#include <dumux/io/format.hh>

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief Symmetrize Dirichlet boundary conditions (recursively) of a coupled problem
 * \note this hard-codes for MultiTypeBlockMatrix<BCRSMatrix<FieldMatrix>> nesting
 * \todo Should be possible to do this generically and recursively with some helper functions
 */
template<class Matrix, class Vector>
void symmetrizeDirichlet(Matrix& A, Vector& b, const Vector& dirichletDofs)
{
    using namespace Dune::Hybrid;
    forEach(std::make_index_sequence<Matrix::N()>{}, [&](const auto i)
    {
        forEach(std::make_index_sequence<Matrix::M()>{}, [&](const auto j)
        {
            auto& level0 = A[i][j];
            auto& level0BRow = b[i];
            const auto& level0BCol = b[j];
            const auto& level0Dirichlet = dirichletDofs[j];
            const bool isDiagonal0 = (i == j);

            const auto rowEnd = level0.end();
            for (auto row = level0.begin(); row != rowEnd; ++row)
            {
                const auto colEnd = row->end();
                for (auto col = row->begin(); col != colEnd; ++col)
                {
                    auto& level1 = level0[row.index()][col.index()];
                    auto& level1BRow = level0BRow[row.index()];
                    const auto& level1BCol = level0BCol[col.index()];
                    const auto& level1Dirichlet = level0Dirichlet[col.index()];
                    const bool isDiagonal1 = ((row.index() == col.index()) && isDiagonal0);

                    for (int jj = 0; jj < level1.M(); ++jj)
                    {
                        // all rows in this column can be eliminated by a Dirichlet row
                        const auto& level2Dirichlet = level1Dirichlet[jj];
                        if (level2Dirichlet > 0.5)
                        {
                            for (int ii = 0; ii < level1.N(); ++ii)
                            {
                                auto& level2 = level1[ii][jj];
                                auto& level2BRow = level1BRow[ii];
                                const auto& level2BCol = level1BCol[jj];
                                const bool isDiagonal2 = ((ii == jj) && isDiagonal1);

                                if (!isDiagonal2)
                                {
                                    level2BRow -= level2*level2BCol;
                                    level2 = 0.0;
                                }
                            }
                        }
                    }
                }
            }
        });
    });
}

/*!
 * \ingroup Linear
 * \brief Check if a matrix is symmetric
 * \note this hard-codes for MultiTypeBlockMatrix<BCRSMatrix<FieldMatrix>> nesting
 * \todo Should be possible to do this generically and recursively with some helper functions
 */
template<class Matrix>
bool matrixIsSymmetric(const Matrix& A)
{
    bool isSymmetric = true;

    using namespace Dune::Hybrid;
    forEach(std::make_index_sequence<Matrix::N()>{}, [&](const auto i)
    {
        const auto& diagBlockA = A[i][i];
        for (auto rowIt = diagBlockA.begin(), rowEndIt = diagBlockA.end(); rowIt != rowEndIt; ++rowIt)
        {
            for (auto colIt = rowIt->begin(), colEndIt = rowIt->end(); colIt != colEndIt; ++colIt)
            {
                auto& block = *colIt;
                for (int ii = 0; ii < block.N(); ++ii)
                {
                    for (int jj = 0; jj < block.M(); ++jj)
                    {
                        if (Dune::FloatCmp::ne(diagBlockA[rowIt.index()][colIt.index()][ii][jj], diagBlockA[colIt.index()][rowIt.index()][jj][ii], 1e-7))
                        {
                            std::cout << Fmt::format("diagBlockA[{}][{}][{}][{}][{}][{}]: ", i, i, rowIt.index(), colIt.index(), ii, jj)
                                    << diagBlockA[rowIt.index()][colIt.index()][ii][jj] << " <--> " << diagBlockA[colIt.index()][rowIt.index()][jj][ii] << std::endl;
                            isSymmetric = false;
                        }
                    }
                }
            }
        }

        forEach(std::make_index_sequence<Matrix::M()>{}, [&](const auto j)
        {
            if constexpr (i != j && j > i)
            {
                const auto& offDiag = A[i][j];
                const auto& offDiagT = A[j][i];
                for (auto rowIt = offDiag.begin(), rowEndIt = offDiag.end(); rowIt != rowEndIt; ++rowIt)
                    for (auto colIt = rowIt->begin(), colEndIt = rowIt->end(); colIt != colEndIt; ++colIt)
                        if (offDiagT.exists(colIt.index(), rowIt.index()))
                        {
                            auto& block = *colIt;
                            for (int ii = 0; ii < block.N(); ++ii)
                            {
                                for (int jj = 0; jj < block.M(); ++jj)
                                {
                                    if (Dune::FloatCmp::ne(offDiag[rowIt.index()][colIt.index()][ii][jj], offDiagT[colIt.index()][rowIt.index()][jj][ii], 1e-7))
                                    {
                                        std::cout << Fmt::format("offDiag[{}][{}][{}][{}][{}][{}]: ", i, j, rowIt.index(), colIt.index(), ii, jj)
                                                << offDiag[rowIt.index()][colIt.index()][ii][jj]
                                                << " <--> "
                                                << offDiagT[colIt.index()][rowIt.index()][jj][ii]
                                                << Fmt::format(" :offDiag[{}][{}][{}][{}][{}][{}]: ", i, j, colIt.index(), rowIt.index(), jj, ii)
                                                << std::endl;
                                        isSymmetric = false;
                                    }
                                }
                            }
                        }
            }
        });
    });

    return isSymmetric;
}

} // end namespace Dumux

#endif
