// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Assembly
 * \brief Helper function to generate Jacobian pattern for different discretization methods
 */
#ifndef DUMUX_JACOBIAN_PATTERN_HH
#define DUMUX_JACOBIAN_PATTERN_HH

#include <type_traits>
#include <dune/istl/matrixindexset.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup Assembly
 * \brief Helper function to generate Jacobian pattern for cell-centered methods
 */
template<bool isImplicit, class GridGeometry,
         typename std::enable_if_t<( (GridGeometry::discMethod == DiscretizationMethods::cctpfa)
                                     || (GridGeometry::discMethod == DiscretizationMethods::ccmpfa) ), int> = 0>
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

/*!
 * \ingroup Assembly
 * \brief Helper function to generate Jacobian pattern for the staggered method
 */
template<bool isImplicit, class GridGeometry,
         typename std::enable_if_t<( (GridGeometry::discMethod == DiscretizationMethods::staggered) ), int> = 0>
auto getJacobianPattern(const GridGeometry& gridGeometry)
{
    // resize the jacobian and the residual
    const auto numDofs = gridGeometry.numDofs();
    Dune::MatrixIndexSet pattern(numDofs, numDofs);

    const auto& connectivityMap = gridGeometry.connectivityMap();

    auto fvGeometry = localView(gridGeometry);
    // evaluate the actual pattern
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        if(gridGeometry.isCellCenter())
        {
            // the global index of the element at hand
            static constexpr auto cellCenterIdx = GridGeometry::cellCenterIdx();
            const auto ccGlobalI = gridGeometry.elementMapper().index(element);
            pattern.add(ccGlobalI, ccGlobalI);
            for (auto&& ccGlobalJ : connectivityMap(cellCenterIdx, cellCenterIdx, ccGlobalI))
                pattern.add(ccGlobalI, ccGlobalJ);
        }
        else
        {
            static constexpr auto faceIdx = GridGeometry::faceIdx();
            fvGeometry.bindElement(element);

            // loop over sub control faces
            for (auto&& scvf : scvfs(fvGeometry))
            {
                const auto faceGlobalI = scvf.dofIndex();
                pattern.add(faceGlobalI, faceGlobalI);
                for (auto&& faceGlobalJ : connectivityMap(faceIdx, faceIdx, scvf.index()))
                    pattern.add(faceGlobalI, faceGlobalJ);
            }
        }
    }

    return pattern;
}

/*!
 * \ingroup Assembly
 * \brief Helper function to generate Jacobian pattern for finite element scheme
 */
template<class FEBasis>
Dune::MatrixIndexSet getFEJacobianPattern(const FEBasis& feBasis)
{
    const auto numDofs = feBasis.size();

    Dune::MatrixIndexSet pattern;
    pattern.resize(numDofs, numDofs);

    auto localView = feBasis.localView();
    // matrix pattern for implicit Jacobians
    for (const auto& element : elements(feBasis.gridView()))
    {
        localView.bind(element);

        const auto& finiteElement = localView.tree().finiteElement();
        const auto numLocalDofs = finiteElement.localBasis().size();
        for (std::size_t i = 0; i < numLocalDofs; i++)
        {
            const auto dofIdxI = localView.index(i);
            for (std::size_t j = 0; j < numLocalDofs; j++)
            {
                const auto dofIdxJ = localView.index(j);
                pattern.add(dofIdxI, dofIdxJ);
            }
        }
    }

    return pattern;
}

/*!
 * \ingroup Assembly
 * \brief Helper function to generate Jacobian pattern for finite element scheme
 * \note This interface is for compatibility with the other schemes. The pattern
 *       in fem is the same independent of the time discretization scheme.
 */
template<bool isImplicit, class GridGeometry,
         typename std::enable_if_t<(GridGeometry::discMethod == DiscretizationMethods::fem), int> = 0>
Dune::MatrixIndexSet getJacobianPattern(const GridGeometry& gridGeometry)
{ return getFEJacobianPattern(gridGeometry.feBasis()); }

/*!
 * \ingroup Assembly
 * \brief Helper function to generate Jacobian pattern for the face-centered methods
 */
template<bool isImplicit, class GridGeometry,
         typename std::enable_if_t<(GridGeometry::discMethod == DiscretizationMethods::fcstaggered), int> = 0>
Dune::MatrixIndexSet getJacobianPattern(const GridGeometry& gridGeometry)
{
    // resize the jacobian and the residual
    const auto numDofs = gridGeometry.numDofs();
    Dune::MatrixIndexSet pattern(numDofs, numDofs);

    const auto& connectivityMap = gridGeometry.connectivityMap();
    auto fvGeometry = localView(gridGeometry);

    // set the pattern
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        fvGeometry.bind(element);
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto globalI = scv.dofIndex();
            pattern.add(globalI, globalI);

            for (const auto& scvIdxJ : connectivityMap[scv.index()])
            {
                const auto globalJ = fvGeometry.scv(scvIdxJ).dofIndex();
                pattern.add(globalI, globalJ);

                if (gridGeometry.isPeriodic())
                {
                    if (gridGeometry.dofOnPeriodicBoundary(globalI) && globalI != globalJ)
                    {
                        const auto globalIP = gridGeometry.periodicallyMappedDof(globalI);
                        pattern.add(globalIP, globalI);
                        pattern.add(globalI, globalIP);

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
 * \ingroup Assembly
 * \brief Compute the Jacobian matrix sparsity pattern for control volume finite element schemes
 */
template<bool isImplicit, class GridGeometry,
         typename std::enable_if_t<DiscretizationMethods::isCVFE<typename GridGeometry::DiscretizationMethod>, int> = 0>
Dune::MatrixIndexSet getJacobianPattern(const GridGeometry& gridGeometry)
{
    // resize the jacobian and the residual
    const auto numDofs = gridGeometry.numDofs();
    Dune::MatrixIndexSet pattern(numDofs, numDofs);

    // matrix pattern for implicit Jacobians
    if constexpr (isImplicit)
    {
        auto fvGeometry = localView(gridGeometry);
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            fvGeometry.bindElement(element);
            for (const auto& scv : scvs(fvGeometry))
            {
                const auto globalI = scv.dofIndex();
                for (const auto& scvJ : scvs(fvGeometry))
                {
                    const auto globalJ = scvJ.dofIndex();
                    pattern.add(globalI, globalJ);

                    if (gridGeometry.isPeriodic())
                    {
                        if (gridGeometry.dofOnPeriodicBoundary(globalI) && globalI != globalJ)
                        {
                            const auto globalIP = gridGeometry.periodicallyMappedDof(globalI);
                            pattern.add(globalIP, globalI);
                            pattern.add(globalI, globalIP);

                            if (globalI > globalIP)
                                pattern.add(globalIP, globalJ);
                        }
                    }
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


} // namespace Dumux

#endif
