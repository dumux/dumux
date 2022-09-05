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
 * \brief Predefined linear algebra backends for Dune linear solvers.
 */
#ifndef DUMUX_LINEAR_DUNE_LINEAR_ALGEBRA_BACKEND_HH
#define DUMUX_LINEAR_DUNE_LINEAR_ALGEBRA_BACKEND_HH

#include <concepts>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <dumux/experimental/new_assembly/dumux/common/concepts.hh>
#include <dumux/experimental/new_assembly/dumux/common/indexstrategies.hh>
#include <dumux/experimental/new_assembly/dumux/assembly/matrixpattern.hh>

// include helper functions for linear algebra backends
#include <dumux/experimental/new_assembly/dumux/linear/labackend.hh>

// include traits for dune linear systems
#include <dumux/experimental/new_assembly/dumux/linear/dune/systemtraits.hh>

namespace Dumux {

template<typename BT, typename A>
void setSparsityPattern(Dune::BCRSMatrix<BT, A>& m, const MatrixPattern& p)
{
    m.setSize(p.numRows(), p.numCols());
    m.setBuildMode(Dune::BCRSMatrix<BT, A>::random);

    for (std::size_t i = 0; i < p.numRows(); ++i)
        m.setrowsize(i, p.rowSize(i));
    m.endrowsizes();

    for (const auto& index : indices(p))
        m.addindex(index.rowIndex, index.colIndex);
    m.endindices();
}

/*!
 * \ingroup Linear
 * \brief Default linear system for Dune linear solvers.
 * \tparam Scalar The type used for field values
 * \tparam n The number of equations per degree of freedom
 */
template<Concepts::Arithmetic Scalar, int n>
struct DefaultBlockedLinearAlgebraBackend
{
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<Scalar, n, n>>;
    using Vector = Dune::BlockVector<Dune::FieldVector<Scalar, n>>;
    using IndexStrategy = BlockedIndexStrategy<2>;

    static IndexStrategy makeIndexStrategy(std::size_t size)
    { return {}; }

    static void setSparsityPattern(Matrix& m, const MatrixPattern& pattern)
    { Dumux::setSparsityPattern(m, pattern); }

    static void setSize(Vector& v, std::size_t size)
    { v.resize(size); }
};

/*!
 * \ingroup Linear
 * \brief Default linear system for Dune linear solvers.
 * \tparam Scalar The type used for field values
 * \tparam n The number of equations per degree of freedom
 */
template<Concepts::Arithmetic Scalar, int n>
struct DefaultFlatLinearAlgebraBackend
{
    using Matrix = Dune::BCRSMatrix<Scalar>;
    using Vector = Dune::BlockVector<Scalar>;
    using IndexStrategy = FlatIndexStrategy<2>;

    static IndexStrategy makeIndexStrategy(std::size_t size)
    { return IndexStrategy{BlockSizes{size, n}}; }

    static void setSparsityPattern(Matrix& m, const MatrixPattern& pattern)
    { Dumux::setSparsityPattern(m, asFlatPattern<n>(pattern)); }

    static void setSize(Vector& v, std::size_t size)
    { v.resize(size); }
};

} // namespace Dumux

#endif
