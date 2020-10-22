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
 * \brief Trilinos solver backend
 */
#ifndef DUMUX_TRILINOS_SOLVER_BACKEND_HH
#define DUMUX_TRILINOS_SOLVER_BACKEND_HH

#include <dune/istl/solver.hh>

#include <dumux/common/typetraits/matrix.hh>
#include <dumux/linear/solver.hh>

#if HAVE_TRILINOS
#include <Tpetra_Vector_decl.hpp>
#include <Tpetra_CrsMatrix_decl.hpp>
#endif

namespace Dumux {

#if HAVE_TRILINOS
/*!
 * \ingroup Linear
 * \brief Linear solver using Trilinos.
 */
class TrilinosSolverBackend : public LinearSolver
{
public:
    using LinearSolver::LinearSolver;

    template<class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        static_assert(isBCRSMatrix<Matrix>::value, "SuperLU only works with BCRS matrices!");
        using BlockType = typename Matrix::block_type;
        static_assert(BlockType::rows == BlockType::cols, "Matrix block must be quadratic!");
        constexpr auto blockSize = BlockType::rows;

        Tpetra::Vector xt;

        for (auto blockIdx = 0u; blockIdx < x.size(); ++blockIdx)
        {
            std::cout << "block idx " << blockIdx << ": ";

            std::cout << "x = ";
            for (auto i = 0u; i < blockSize; ++i)
                std::cout << x[blockIdx][i] << ", ";

            std::cout << "b = ";
            for (auto i = 0u; i < blockSize; ++i)
                std::cout << b[blockIdx][i] << ", ";

            std::cout << std::endl;
        }

        // Tpetra::CrsMatrix At;

        for (auto rowIt = A.begin(); rowIt != A.end(); ++rowIt)
        {
            std::cout << "block row idx " << rowIt.index() << ":" << std::endl;

            for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt)
            {
                std::cout << "  block col idx " << colIt.index() << ": block entry = { ";

                for (auto i = 0u; i < blockSize; ++i)
                {
                    std::cout << "{ ";

                    for (auto j = 0u; j < blockSize; ++j)
                    {
                        std::cout << (*colIt)[i][j] << ", ";
                    }

                    std::cout << " }, ";
                }

                std::cout << " }" << std::endl;
            }
        }

        return result_.converged;
    }

    std::string name() const
    {
        return "Trilinos solver";
    }

    const Dune::InverseOperatorResult& result() const
    {
        return result_;
    }

private:
    Dune::InverseOperatorResult result_;
};
#endif // HAVE_TRILINOS

} // end namespace Dumux

#endif
