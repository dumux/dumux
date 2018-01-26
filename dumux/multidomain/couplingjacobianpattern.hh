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
 * \brief Helper function to generate Jacobian pattern for multidomain models
 */
#ifndef DUMUX_MUTLIDOMAIN_COUPLING_JACOBIAN_PATTERN_HH
#define DUMUX_MUTLIDOMAIN_COUPLING_JACOBIAN_PATTERN_HH

#include <type_traits>
#include <dune/istl/matrixindexset.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux {

/*!
 * \ingroup Assembly
 * \brief Helper function to generate coupling Jacobian pattern (off-diagonal blocks)
 *        for cell-centered schemes
 */
template<bool isImplicit, class CouplingManager, class GridGeometry0, class GridGeometry1, std::size_t i, std::size_t j,
         typename std::enable_if_t<( (GridGeometry0::discretizationMethod == DiscretizationMethods::CCTpfa)
                                     || (GridGeometry0::discretizationMethod == DiscretizationMethods::CCMpfa) ), int> = 0>
Dune::MatrixIndexSet getCouplingJacobianPattern(const CouplingManager& couplingManager,
                                                Dune::index_constant<i> domainI,
                                                const GridGeometry0& gridGeometry0,
                                                Dune::index_constant<j> domainJ,
                                                const GridGeometry1& gridGeometry1)
{
    const auto numDofs0 = gridGeometry0.numDofs();
    const auto numDofs1 = gridGeometry1.numDofs();
    Dune::MatrixIndexSet pattern;
    pattern.resize(numDofs0, numDofs1);

    // matrix pattern for implicit Jacobians
    if (isImplicit)
    {
        for (const auto& element0 : elements(gridGeometry0.gridView()))
        {
            const auto& stencil = couplingManager.couplingStencil(element0, domainI, domainJ);
            const auto globalI = gridGeometry0.elementMapper().index(element0);
            for (const auto globalJ : stencil)
                pattern.add(globalI, globalJ);
        }
    }

    // matrix pattern for explicit Jacobians
    // -> diagonal matrix, so coupling block is empty
    // just return the empty pattern

    return pattern;
}

/*!
 * \ingroup Assembly
 * \brief Helper function to generate coupling Jacobian pattern (off-diagonal blocks)
 *        for the box scheme
 */
template<bool isImplicit, class CouplingManager, class GridGeometry0, class GridGeometry1, std::size_t i, std::size_t j,
         typename std::enable_if_t<(GridGeometry0::discretizationMethod == DiscretizationMethods::Box), int> = 0>
Dune::MatrixIndexSet getCouplingJacobianPattern(const CouplingManager& couplingManager,
                                                Dune::index_constant<i> domainI,
                                                const GridGeometry0& gridGeometry0,
                                                Dune::index_constant<j> domainJ,
                                                const GridGeometry1& gridGeometry1)
{
    const auto numDofs0 = gridGeometry0.numDofs();
    const auto numDofs1 = gridGeometry1.numDofs();
    Dune::MatrixIndexSet pattern;
    pattern.resize(numDofs0, numDofs1);

    // matrix pattern for implicit Jacobians
    if (isImplicit)
    {
        static constexpr int dim = std::decay_t<decltype(gridGeometry0.gridView())>::dimension;
        for (const auto& element0 : elements(gridGeometry0.gridView()))
        {
            const auto& stencil = couplingManager.couplingStencil(element0, domainI, domainJ);
            for (std::size_t vIdxLocal = 0; vIdxLocal < element0.subEntities(dim); ++vIdxLocal)
            {
                const auto globalI = gridGeometry0.vertexMapper().subIndex(element0, vIdxLocal, dim);
                for (const auto globalJ : stencil)
                    pattern.add(globalI, globalJ);
            }
        }
    }

    // matrix pattern for explicit Jacobians
    // -> diagonal matrix, so coupling block is empty
    // just return the empty pattern

    return pattern;
}

} // namespace Dumux

#endif
