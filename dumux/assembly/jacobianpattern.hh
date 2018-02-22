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
 * \ingroup Assembly
 * \brief Helper function to generate Jacobian pattern for different discretization methods
 */
#ifndef DUMUX_JACOBIAN_PATTERN_HH
#define DUMUX_JACOBIAN_PATTERN_HH

#include <type_traits>
#include <dune/istl/matrixindexset.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux {

/*!
 * \ingroup Assembly
 * \brief Helper function to generate Jacobian pattern for the box method
 */
template<bool isImplicit, class GridGeometry,
         typename std::enable_if_t<(GridGeometry::discretizationMethod == DiscretizationMethod::box), int> = 0>
Dune::MatrixIndexSet getJacobianPattern(const GridGeometry& gridGeometry)
{
    const auto numDofs = gridGeometry.numDofs();
    Dune::MatrixIndexSet pattern;
    pattern.resize(numDofs, numDofs);

    // matrix pattern for implicit Jacobians
    if (isImplicit)
    {
        static constexpr int dim = std::decay_t<decltype(gridGeometry.gridView())>::dimension;
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            for (unsigned int vIdx = 0; vIdx < element.subEntities(dim); ++vIdx)
            {
                const auto globalI = gridGeometry.vertexMapper().subIndex(element, vIdx, dim);
                for (unsigned int vIdx2 = vIdx; vIdx2 < element.subEntities(dim); ++vIdx2)
                {
                    const auto globalJ = gridGeometry.vertexMapper().subIndex(element, vIdx2, dim);
                    pattern.add(globalI, globalJ);
                    pattern.add(globalJ, globalI);
                }
            }
        }
    }

    // matrix pattern for explicit Jacobians -> diagonal matrix
    else
    {
        for (unsigned int globalI = 0; globalI < numDofs; ++globalI)
            pattern.add(globalI, globalI);
    }

    return pattern;
}

/*!
 * \ingroup Assembly
 * \brief Helper function to generate Jacobian pattern for cell-centered methods
 */
template<bool isImplicit, class GridGeometry,
         typename std::enable_if_t<( (GridGeometry::discretizationMethod == DiscretizationMethod::cctpfa)
                                     || (GridGeometry::discretizationMethod == DiscretizationMethod::ccmpfa) ), int> = 0>
Dune::MatrixIndexSet getJacobianPattern(const GridGeometry& gridGeometry)
{
    const auto numDofs = gridGeometry.numDofs();
    Dune::MatrixIndexSet pattern;
    pattern.resize(numDofs, numDofs);

    // matrix pattern for implicit Jacobians
    if (isImplicit)
    {
        for (unsigned int globalI = 0; globalI < numDofs; ++globalI)
        {
            pattern.add(globalI, globalI);
            for (const auto& dataJ : gridGeometry.connectivityMap()[globalI])
                pattern.add(dataJ.globalJ, globalI);
        }
    }

    // matrix pattern for explicit Jacobians -> diagonal matrix
    else
    {
        for (unsigned int globalI = 0; globalI < numDofs; ++globalI)
            pattern.add(globalI, globalI);
    }

    return pattern;
}

} // namespace Dumux

#endif
