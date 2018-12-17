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
 * \ingroup MultiDomain
 * \brief Helper function to generate Jacobian pattern for multi domain models
 */
#ifndef DUMUX_MUTLIDOMAIN_COUPLING_JACOBIAN_PATTERN_HH
#define DUMUX_MUTLIDOMAIN_COUPLING_JACOBIAN_PATTERN_HH

#include <type_traits>
#include <dune/common/indices.hh>
#include <dune/istl/matrixindexset.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \brief Helper function to generate coupling Jacobian pattern (off-diagonal blocks)
 *        for cell-centered schemes
 */
template<bool isImplicit, class CouplingManager, class GridGeometryI, class GridGeometryJ, std::size_t i, std::size_t j,
         typename std::enable_if_t<( (GridGeometryI::discMethod == DiscretizationMethod::cctpfa)
                                     || (GridGeometryI::discMethod == DiscretizationMethod::ccmpfa) ), int> = 0>
Dune::MatrixIndexSet getCouplingJacobianPattern(const CouplingManager& couplingManager,
                                                Dune::index_constant<i> domainI,
                                                const GridGeometryI& gridGeometryI,
                                                Dune::index_constant<j> domainJ,
                                                const GridGeometryJ& gridGeometryJ)
{
    const auto numDofsI = gridGeometryI.numDofs();
    const auto numDofsJ = gridGeometryJ.numDofs();
    Dune::MatrixIndexSet pattern;
    pattern.resize(numDofsI, numDofsJ);

    // matrix pattern for implicit Jacobians
    if (isImplicit)
    {
        for (const auto& elementI : elements(gridGeometryI.gridView()))
        {
            const auto& stencil = couplingManager.couplingStencil(domainI, elementI, domainJ);
            const auto globalI = gridGeometryI.elementMapper().index(elementI);
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
 * \ingroup MultiDomain
 * \brief Helper function to generate coupling Jacobian pattern (off-diagonal blocks)
 *        for the box scheme
 */
template<bool isImplicit, class CouplingManager, class GridGeometryI, class GridGeometryJ, std::size_t i, std::size_t j,
         typename std::enable_if_t<(GridGeometryI::discMethod == DiscretizationMethod::box), int> = 0>
Dune::MatrixIndexSet getCouplingJacobianPattern(const CouplingManager& couplingManager,
                                                Dune::index_constant<i> domainI,
                                                const GridGeometryI& gridGeometryI,
                                                Dune::index_constant<j> domainJ,
                                                const GridGeometryJ& gridGeometryJ)
{
    const auto numDofsI = gridGeometryI.numDofs();
    const auto numDofsJ = gridGeometryJ.numDofs();
    Dune::MatrixIndexSet pattern;
    pattern.resize(numDofsI, numDofsJ);

    // matrix pattern for implicit Jacobians
    if (isImplicit)
    {
        static constexpr int dim = std::decay_t<decltype(gridGeometryI.gridView())>::dimension;
        for (const auto& elementI : elements(gridGeometryI.gridView()))
        {
            const auto& stencil = couplingManager.couplingStencil(domainI, elementI, domainJ);
            for (std::size_t vIdxLocal = 0; vIdxLocal < elementI.subEntities(dim); ++vIdxLocal)
            {
                const auto globalI = gridGeometryI.vertexMapper().subIndex(elementI, vIdxLocal, dim);
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

/*!
 * \ingroup MultiDomain
 * \brief Helper function to generate coupling Jacobian pattern (off-diagonal blocks)
 *        for the staggered scheme (degrees of freedom on cell centers)
 */
template<bool isImplicit, class CouplingManager, class GridGeometryI, class GridGeometryJ, std::size_t i, std::size_t j,
         typename std::enable_if_t<(GridGeometryI::discMethod == DiscretizationMethod::staggered &&
                                    GridGeometryI::isCellCenter()), int> = 0>
Dune::MatrixIndexSet getCouplingJacobianPattern(const CouplingManager& couplingManager,
                                                Dune::index_constant<i> domainI,
                                                const GridGeometryI& gridGeometryI,
                                                Dune::index_constant<j> domainJ,
                                                const GridGeometryJ& gridGeometryJ)
{
    Dune::MatrixIndexSet pattern(gridGeometryI.numDofs(), gridGeometryJ.numDofs());

    for (const auto& elementI : elements(gridGeometryI.gridView()))
    {
        const auto ccGlobalI = gridGeometryI.elementMapper().index(elementI);
        for (auto&& faceGlobalJ : couplingManager.couplingStencil(domainI, elementI, domainJ))
                pattern.add(ccGlobalI, faceGlobalJ);
    }

    return pattern;
}

/*!
 * \ingroup MultiDomain
 * \brief Helper function to generate coupling Jacobian pattern (off-diagonal blocks)
 *        for the staggered scheme (degrees of freedom on faces)
 */
template<bool isImplicit, class CouplingManager, class GridGeometryI, class GridGeometryJ, std::size_t i, std::size_t j,
         typename std::enable_if_t<(GridGeometryI::discMethod == DiscretizationMethod::staggered &&
                                    GridGeometryI::isFace()), int> = 0>
Dune::MatrixIndexSet getCouplingJacobianPattern(const CouplingManager& couplingManager,
                                                Dune::index_constant<i> domainI,
                                                const GridGeometryI& gridGeometryI,
                                                Dune::index_constant<j> domainJ,
                                                const GridGeometryJ& gridGeometryJ)
{
    Dune::MatrixIndexSet pattern(gridGeometryI.numDofs(), gridGeometryJ.numDofs());

    for (const auto& elementI : elements(gridGeometryI.gridView()))
    {
        auto fvGeometry = localView(gridGeometryI);
        fvGeometry.bindElement(elementI);

        // loop over sub control faces
        for (auto&& scvf : scvfs(fvGeometry))
        {
            const auto faceGlobalI = scvf.dofIndex();
            for (auto&& globalJ : couplingManager.couplingStencil(domainI, scvf, domainJ))
                pattern.add(faceGlobalI, globalJ);
        }
    }

    return pattern;
}

} // end namespace Dumux

#endif
