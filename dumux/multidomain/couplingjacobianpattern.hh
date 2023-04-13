// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
         typename std::enable_if_t<( (GridGeometryI::discMethod == DiscretizationMethods::cctpfa)
                                     || (GridGeometryI::discMethod == DiscretizationMethods::ccmpfa) ), int> = 0>
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
 *        for the staggered scheme (degrees of freedom on cell centers)
 */
template<bool isImplicit, class CouplingManager, class GridGeometryI, class GridGeometryJ, std::size_t i, std::size_t j,
         typename std::enable_if_t<(GridGeometryI::discMethod == DiscretizationMethods::staggered &&
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
         typename std::enable_if_t<(GridGeometryI::discMethod == DiscretizationMethods::staggered &&
                                    GridGeometryI::isFace()), int> = 0>
Dune::MatrixIndexSet getCouplingJacobianPattern(const CouplingManager& couplingManager,
                                                Dune::index_constant<i> domainI,
                                                const GridGeometryI& gridGeometryI,
                                                Dune::index_constant<j> domainJ,
                                                const GridGeometryJ& gridGeometryJ)
{
    Dune::MatrixIndexSet pattern(gridGeometryI.numDofs(), gridGeometryJ.numDofs());

    auto fvGeometry = localView(gridGeometryI);
    for (const auto& elementI : elements(gridGeometryI.gridView()))
    {
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

/*!
 * \ingroup MultiDomain
 * \brief Helper function to generate coupling Jacobian pattern (off-diagonal blocks)
 *        for the staggered scheme (degrees of freedom on cell centers)
 */
template<bool isImplicit, class CouplingManager, class GridGeometryI, class GridGeometryJ, std::size_t i, std::size_t j,
         typename std::enable_if_t<(GridGeometryI::discMethod == DiscretizationMethods::fcstaggered), int> = 0>
Dune::MatrixIndexSet getCouplingJacobianPattern(const CouplingManager& couplingManager,
                                                Dune::index_constant<i> domainI,
                                                const GridGeometryI& gridGeometryI,
                                                Dune::index_constant<j> domainJ,
                                                const GridGeometryJ& gridGeometryJ)
{
    Dune::MatrixIndexSet pattern(gridGeometryI.numDofs(), gridGeometryJ.numDofs());

    auto fvGeometry = localView(gridGeometryI);
    for (const auto& elementI : elements(gridGeometryI.gridView()))
    {
        fvGeometry.bindElement(elementI);
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto globalI = scv.dofIndex();
            const auto& stencil = couplingManager.couplingStencil(domainI, elementI, scv, domainJ);
            for (const auto globalJ : stencil)
            {
                assert(globalJ < gridGeometryJ.numDofs());
                pattern.add(globalI, globalJ);

                if (gridGeometryI.isPeriodic())
                {
                    if (gridGeometryI.dofOnPeriodicBoundary(globalI))
                    {
                        const auto globalIP = gridGeometryI.periodicallyMappedDof(globalI);

                        if (globalI > globalIP)
                            pattern.add(globalIP, globalJ);
                    }
                }
            }
        }
    }

    return pattern;
}

/*!
 * \ingroup MultiDomain
 * \brief Helper function to generate coupling Jacobian pattern (off-diagonal blocks)
 *        for control-volume finite element schemes (CVFE)
 */
template<bool isImplicit, class CouplingManager, class GridGeometryI, class GridGeometryJ, std::size_t i, std::size_t j,
         typename std::enable_if_t<DiscretizationMethods::isCVFE<typename GridGeometryI::DiscretizationMethod>, int> = 0>
Dune::MatrixIndexSet getCouplingJacobianPattern(const CouplingManager& couplingManager,
                                                Dune::index_constant<i> domainI,
                                                const GridGeometryI& gridGeometryI,
                                                Dune::index_constant<j> domainJ,
                                                const GridGeometryJ& gridGeometryJ)
{
    Dune::MatrixIndexSet pattern;

    // matrix pattern for implicit Jacobians
    if constexpr (isImplicit)
    {
        pattern.resize(gridGeometryI.numDofs(),  gridGeometryJ.numDofs());
        auto fvGeometry = localView(gridGeometryI);
        for (const auto& elementI : elements(gridGeometryI.gridView()))
        {
            fvGeometry.bindElement(elementI);
            const auto& stencil = couplingManager.couplingStencil(domainI, elementI, domainJ);
            for (const auto& scv : scvs(fvGeometry))
            {
                for (const auto globalJ : stencil)
                    pattern.add(scv.dofIndex(), globalJ);

            }
        }
    }

    // matrix pattern for explicit Jacobians
    // -> diagonal matrix, so coupling block is empty
    // just return the empty pattern

    return pattern;
}

} // end namespace Dumux

#endif
