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
 * \brief Concepts and helper functions for linear algebra backends.
 */
#ifndef DUMUX_LINEAR_LA_BACKEND_HH
#define DUMUX_LINEAR_LA_BACKEND_HH

#include <memory>
#include <concepts>

#include <dumux/experimental/new_assembly/dumux/linear/system.hh>
#include <dumux/experimental/new_assembly/dumux/linear/concepts.hh>

#include <dumux/experimental/new_assembly/dumux/assembly/matrixpattern.hh>
#include <dumux/experimental/new_assembly/dumux/assembly/jacobianpattern.hh>

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief Helper function to make a zero initialized matrix
 *        from a linear algebra backend and a given sparsity pattern.
 */
template<Concepts::LinearAlgebraBackend LAB> requires(
    std::default_initializable<typename LAB::Matrix>)
typename LAB::Matrix makeMatrix(const LAB& backend,
                                const MatrixPattern& pattern)
{
    typename LAB::Matrix matrix;
    backend.setSparsityPattern(matrix, pattern);
    LinearSystem::fill(matrix, 0.0);
    return matrix;
}


/*!
 * \ingroup Linear
 * \brief Helper function to make a zero initialized matrix
 *        from a linear algebra backend and a given grid geometry.
 */
template<Concepts::LinearAlgebraBackend LAB, typename GridGeometry>
typename LAB::Matrix makeMatrix(const LAB& backend,
                                const GridGeometry& gridGeometry)
{ return makeMatrix(backend, getJacobianPattern(gridGeometry)); }


/*!
 * \ingroup Linear
 * \brief Helper function to make a zero-initialized vector
 *        from a linear algebra backend and given size.
 */
template<Concepts::LinearAlgebraBackend LAB> requires(
    std::default_initializable<typename LAB::Vector>)
typename LAB::Vector makeVector(const LAB& backend, std::size_t size)
{
    typename LAB::Vector vector;
    backend.setSize(vector, size);
    LinearSystem::fill(vector, 0.0);
    return vector;
}


/*!
 * \ingroup Linear
 * \brief Helper function to make a zero-initialized vector
 *        from a linear algebra backend and given grid geometry.
 */
template<Concepts::LinearAlgebraBackend LAB, typename GridGeometry>
auto makeVector(const LAB& backend, const GridGeometry& gridGeometry)
{ return makeVector(backend, gridGeometry.numDofs()); }

} // namespace Dumux

#endif
