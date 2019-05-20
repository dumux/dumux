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
 * \brief Helper function to generate Jacobian pattern for the box method
 */
template<bool isImplicit, class GridGeometry,
         typename std::enable_if_t<(GridGeometry::discMethod == DiscretizationMethod::box), int> = 0>
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
                for (unsigned int vIdx2 = 0; vIdx2 < element.subEntities(dim); ++vIdx2)
                {
                    const auto globalJ = gridGeometry.vertexMapper().subIndex(element, vIdx2, dim);
                    pattern.add(globalI, globalJ);

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
         typename std::enable_if_t<( (GridGeometry::discMethod == DiscretizationMethod::cctpfa)
                                     || (GridGeometry::discMethod == DiscretizationMethod::ccmpfa) ), int> = 0>
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
         typename std::enable_if_t<( (GridGeometry::discMethod == DiscretizationMethod::staggered) ), int> = 0>
auto getJacobianPattern(const GridGeometry& gridGeometry)
{
    // resize the jacobian and the residual
    const auto numDofs = gridGeometry.numDofs();
    Dune::MatrixIndexSet pattern(numDofs, numDofs);

    const auto& connectivityMap = gridGeometry.connectivityMap();

    // evaluate the acutal pattern
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
            auto fvGeometry = localView(gridGeometry);
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

    // matrix pattern for implicit Jacobians
    for (const auto& element : elements(feBasis.gridView()))
    {
        auto localView = feBasis.localView();
        localView.bind(element);

        const auto& finiteElement = localView.tree().finiteElement();
        const auto numLocalDofs = finiteElement.localBasis().size();
        for (size_t i = 0; i < numLocalDofs; i++)
        {
            const auto dofIdxI = localView.index(i);
            for (size_t j = 0; j < numLocalDofs; j++)
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
         typename std::enable_if_t<(GridGeometry::discMethod == DiscretizationMethod::fem), int> = 0>
Dune::MatrixIndexSet getJacobianPattern(const GridGeometry& gridGeometry)
{ return getFEJacobianPattern(gridGeometry.feBasis()); }

} // namespace Dumux

#endif
