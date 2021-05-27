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
 * \brief Coloring schemes for shared-memory-parallel assembly
 */
#ifndef DUMUX_ASSEMBLY_COLORING_HH
#define DUMUX_ASSEMBLY_COLORING_HH

#include <vector>
#include <deque>
#include <iostream>
#include <tuple>

#include <dune/common/timer.hh>
#include <dune/common/exceptions.hh>

#include <dumux/io/format.hh>
#include <dumux/discretization/method.hh>
#include <dumux/assembly/coloring.hh>

namespace Dumux::Detail {

//! Compute a map from dof indices to element indices (helper data for coloring algorithm)
template <std::size_t i, std::size_t j, class CouplingManager>
std::vector<std::vector<std::size_t>>
computeDofToElementMap(Dune::index_constant<i>, Dune::index_constant<j>,
                       const CouplingManager& couplingManager)
{
    std::vector<std::vector<std::size_t>> dofToElements;

    if constexpr (GridGeometry::discMethod == DiscretizationMethod::cctpfa)
    {
        dofToElements.resize(gg.gridView().size(0));
        const auto& eMapper = gg.elementMapper();
        for (const auto& element : elements(gg.gridView()))
        {
            const auto eIdx = eMapper.index(element);
            for (const auto& intersection : intersections(gg.gridView(), element))
                if (intersection.neighbor())
                    dofToElements[eMapper.index(intersection.outside())].push_back(eIdx);
        }
    }

    else if constexpr (GridGeometry::discMethod == DiscretizationMethod::box)
    {
        static constexpr int dim = GridGeometry::GridView::dimension;
        dofToElements.resize(gg.gridView().size(dim));
        const auto& vMapper = gg.vertexMapper();
        for (const auto& element : elements(gg.gridView()))
        {
            const auto eIdx = gg.elementMapper().index(element);
            for (int i = 0; i < element.subEntities(dim); i++)
                dofToElements[vMapper.subIndex(element, i, dim)].push_back(eIdx);
        }
    }

    else
        DUNE_THROW(Dune::NotImplemented,
            "Missing coloring scheme implementation for this discretization method");

    return dofToElements;
}

template <std::size_t i, class CouplingManager>
std::vector<std::vector<std::size_t>>
computeDofToElementMap(Dune::index_constant<i> domainI, Dune::index_constant<i> domainJ,
                       const CouplingManager& couplingManager)
{
    // get regular (single domain) map
    const auto& gg = couplingManager.problem(domainI).gridGeometry();
    auto dofToElement = Detail::computeDofToElementMap(gg);


}

/*!
 * \brief Compute the colors of neighboring nodes in the dependency graph
 *
 * Neighboring nodes are those elements that manipulate the
 * same data structures (e.g. system matrix) in the same places
 *
 * \param gridGeometry the grid geometry
 * \param element the element we want to color
 * \param colors a vector of current colors for each element (not assigned: -1)
 * \param dofToElement a map from dof indices to element indices
 * \param neighborColors a vector to add the colors of neighbor nodes to
 */
template<class GridGeometry, class DofToElementMap>
void addNeighborColors(const GridGeometry& gg,
                       const typename GridGeometry::LocalView::Element& element,
                       const std::vector<int>& colors,
                       const DofToElementMap& dofToElement,
                       std::vector<int>& neighborColors)
{
    if constexpr (GridGeometry::discMethod == DiscretizationMethod::cctpfa)
    {
        // we modify neighbor elements during the assembly
        // check who else modifies these neighbor elements
        const auto& eMapper = gg.elementMapper();
        for (const auto& intersection : intersections(gg.gridView(), element))
            if (intersection.neighbor())
                for (auto eIdx : dofToElement[eMapper.index(intersection.outside())])
                    neighborColors.push_back(colors[eIdx]);
    }

    else if constexpr (GridGeometry::discMethod == DiscretizationMethod::box)
    {
        // we modify the vertex dofs of our element during the assembly
        // check who else modifies these vertex dofs
        const auto& vMapper = gg.vertexMapper();
        static constexpr int dim = GridGeometry::GridView::dimension;
        for (int i = 0; i < element.subEntities(dim); i++)
            for (auto eIdx : dofToElement[vMapper.subIndex(element, i, dim)])
                neighborColors.push_back(colors[eIdx]);
    }

    else
        DUNE_THROW(Dune::NotImplemented,
            "Missing coloring scheme implementation for this discretization method");
}

} // end namespace Dumux::Detail

#endif
