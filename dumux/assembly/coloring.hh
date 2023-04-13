// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
#include <algorithm>

#include <dune/common/timer.hh>
#include <dune/common/exceptions.hh>

#include <dumux/io/format.hh>
#include <dumux/discretization/method.hh>

#ifndef DOXYGEN // hide from doxygen
namespace Dumux::Detail {

//! Compute a map from dof indices to element indices (helper data for coloring algorithm)
template <class GridGeometry>
std::vector<std::vector<std::size_t>>
computeConnectedElements(const GridGeometry& gg)
{
    std::vector<std::vector<std::size_t>> connectedElements;

    if constexpr (GridGeometry::discMethod == DiscretizationMethods::cctpfa)
    {
        connectedElements.resize(gg.gridView().size(0));
        const auto& eMapper = gg.elementMapper();
        for (const auto& element : elements(gg.gridView()))
        {
            const auto eIdx = eMapper.index(element);
            for (const auto& intersection : intersections(gg.gridView(), element))
                if (intersection.neighbor())
                    connectedElements[eMapper.index(intersection.outside())].push_back(eIdx);
        }
    }

    else if constexpr (GridGeometry::discMethod == DiscretizationMethods::box
                       || GridGeometry::discMethod == DiscretizationMethods::pq1bubble)
    {
        static constexpr int dim = GridGeometry::GridView::dimension;
        connectedElements.resize(gg.gridView().size(dim));
        const auto& vMapper = gg.vertexMapper();
        for (const auto& element : elements(gg.gridView()))
        {
            const auto eIdx = gg.elementMapper().index(element);
            for (int i = 0; i < element.subEntities(dim); i++)
                connectedElements[vMapper.subIndex(element, i, dim)].push_back(eIdx);
        }
    }

    else if constexpr (
        GridGeometry::discMethod == DiscretizationMethods::fcstaggered
        || GridGeometry::discMethod == DiscretizationMethods::ccmpfa
    )
    {
        // for MPFA-O schemes the assembly of each element residual touches all vertex neighbors
        // for face-centered staggered it is all codim-2 neighbors (vertex neighbors in 2D, edge neighbors in 3D)
        // but we use vertex neighbors also in 3D for simplicity
        std::vector<std::vector<std::size_t>> vToElements;
        static constexpr int dim = GridGeometry::GridView::dimension;
        vToElements.resize(gg.gridView().size(dim));
        const auto& vMapper = gg.vertexMapper();
        for (const auto& element : elements(gg.gridView()))
        {
            const auto eIdx = gg.elementMapper().index(element);
            for (int i = 0; i < element.subEntities(dim); i++)
                vToElements[vMapper.subIndex(element, i, dim)].push_back(eIdx);
        }

        connectedElements.resize(gg.gridView().size(0));
        for (const auto& element : elements(gg.gridView()))
        {
            const auto eIdx = gg.elementMapper().index(element);
            for (int i = 0; i < element.subEntities(dim); i++)
            {
                const auto& e = vToElements[vMapper.subIndex(element, i, dim)];
                connectedElements[eIdx].insert(connectedElements[eIdx].end(), e.begin(), e.end());
            }

            // make unique
            std::sort(connectedElements[eIdx].begin(), connectedElements[eIdx].end());
            connectedElements[eIdx].erase(
                std::unique(connectedElements[eIdx].begin(), connectedElements[eIdx].end()),
                connectedElements[eIdx].end()
            );
        }
    }

    // nothing has to be precomputed here as only immediate face neighbors are connected
    else if constexpr (GridGeometry::discMethod == DiscretizationMethods::fcdiamond)
        return connectedElements;

    else
        DUNE_THROW(Dune::NotImplemented,
            "Missing coloring scheme implementation for this discretization method"
        );

    return connectedElements;
}

/*!
 * \brief Compute the colors of neighboring nodes in the dependency graph
 *
 * Neighboring nodes are those elements that manipulate the
 * same data structures (e.g. system matrix, volvars, flux cache) in the same places
 *
 * \param gridGeometry the grid geometry
 * \param element the element we want to color
 * \param colors a vector of current colors for each element (not assigned: -1)
 * \param connectedElements a map from implementation-defined indices to element indices
 * \param neighborColors a vector to add the colors of neighbor nodes to
 */
template<class GridGeometry, class ConnectedElements>
void addNeighborColors(const GridGeometry& gg,
                       const typename GridGeometry::LocalView::Element& element,
                       const std::vector<int>& colors,
                       const ConnectedElements& connectedElements,
                       std::vector<int>& neighborColors)
{
    if constexpr (
        GridGeometry::discMethod == DiscretizationMethods::cctpfa
        || GridGeometry::discMethod == DiscretizationMethods::ccmpfa
        || GridGeometry::discMethod == DiscretizationMethods::fcstaggered
    )
    {
        // we modify neighbor elements during the assembly
        // check who else modifies these neighbor elements
        const auto& eMapper = gg.elementMapper();
        for (const auto& intersection : intersections(gg.gridView(), element))
        {
            if (intersection.neighbor())
            {
                // direct face neighbors
                const auto nIdx = eMapper.index(intersection.outside());
                neighborColors.push_back(colors[nIdx]);

                // neighbor-neighbors
                for (const auto nnIdx : connectedElements[eMapper.index(intersection.outside())])
                    neighborColors.push_back(colors[nnIdx]);
            }
        }
    }

    else if constexpr (GridGeometry::discMethod == DiscretizationMethods::box
                       || GridGeometry::discMethod == DiscretizationMethods::pq1bubble)
    {
        // we modify the vertex dofs of our element during the assembly
        // check who else modifies these vertex dofs
        const auto& vMapper = gg.vertexMapper();
        static constexpr int dim = GridGeometry::GridView::dimension;
        // direct vertex neighbors
        for (int i = 0; i < element.subEntities(dim); i++)
            for (auto eIdx : connectedElements[vMapper.subIndex(element, i, dim)])
                neighborColors.push_back(colors[eIdx]);
    }

    else if constexpr (GridGeometry::discMethod == DiscretizationMethods::fcdiamond)
    {
        // we modify neighbor faces during the assembly
        // check who else modifies these neighbor elements
        const auto& eMapper = gg.elementMapper();
        for (const auto& intersection : intersections(gg.gridView(), element))
            if (intersection.neighbor())
                neighborColors.push_back(colors[eMapper.index(intersection.outside())]);
    }

    else
        DUNE_THROW(Dune::NotImplemented,
            "Missing coloring scheme implementation for this discretization method"
        );
}

/*!
 * \brief Find the smallest color (integer >= 0) _not_ present in the given list of colors
 * \param colors list of colors which are already taken
 * \param notAssigned container to store which colors are not yet taken (is resized as required)
 */
int smallestAvailableColor(const std::vector<int>& colors,
                           std::vector<bool>& colorUsed)
{
    const int numColors = colors.size();
    colorUsed.assign(numColors, false);

    // The worst case for e.g. numColors=3 is colors={0, 1, 2}
    // in which case we return 3 as smallest available color
    // That means, we only track candidates in the (half-open) interval [0, numColors)
    // Mark candidate colors which are present in colors
    for (int i = 0; i < numColors; i++)
        if (colors[i] >= 0 && colors[i] < numColors)
            colorUsed[colors[i]] = true;

    // return smallest color not in colors
    for (int i = 0; i < numColors; i++)
        if (!colorUsed[i])
            return i;

    return numColors;
}

} // end namespace Dumux::Detail
#endif // DOXYGEN

namespace Dumux {

/*!
 * \brief Compute iterable lists of element seeds partitioned by color
 *
 * Splits up the elements of a grid view into partitions such that
 * all elements in one partition do not modify global data structures
 * at the same place during assembly. This is used to allow for
 * lock-free thread-parallel (shared memory) assembly routines.
 *
 * Implements a simply greedy graph coloring algorithm:
 * For each node (element), assign the smallest available color
 * not used by any of the neighboring nodes (element with conflicting memory access)
 * The greedy algorithm doesn't necessarily return the smallest
 * possible number of colors (that's a hard problem) but is fast
 *
 * Returns a struct with access to the colors of each element (member colors)
 * and vector of element seed sets of the same color (member sets)
 *
 * \param gg the grid geometry
 * \param verbosity the verbosity level
 */
template<class GridGeometry>
auto computeColoring(const GridGeometry& gg, int verbosity = 1)
{
    Dune::Timer timer;

    using ElementSeed = typename GridGeometry::GridView::Grid::template Codim<0>::EntitySeed;
    struct Coloring
    {
        using Sets = std::deque<std::vector<ElementSeed>>;
        using Colors = std::vector<int>;

        Coloring(std::size_t size) : sets{}, colors(size, -1) {}

        Sets sets;
        Colors colors;
    };

    Coloring coloring(gg.gridView().size(0));

    // pre-reserve some memory for helper arrays to avoid reallocation
    std::vector<int> neighborColors; neighborColors.reserve(30);
    std::vector<bool> colorUsed; colorUsed.reserve(30);

    // dof to element map to speed up neighbor search
    const auto connectedElements = Detail::computeConnectedElements(gg);

    for (const auto& element : elements(gg.gridView()))
    {
        // compute neighbor colors based on discretization-dependent stencil
        neighborColors.clear();
        Detail::addNeighborColors(gg, element, coloring.colors, connectedElements, neighborColors);

        // find smallest color (positive integer) not in neighborColors
        const auto color = Detail::smallestAvailableColor(neighborColors, colorUsed);

        // assign color to element
        coloring.colors[gg.elementMapper().index(element)] = color;

        // add element to the set of elements with the same color
        if (color < coloring.sets.size())
            coloring.sets[color].push_back(element.seed());
        else
            coloring.sets.push_back(std::vector<ElementSeed>{ element.seed() });
    }

    if (verbosity > 0)
        std::cout << Fmt::format("Colored {} elements with {} colors in {} seconds.\n",
                                 gg.gridView().size(0), coloring.sets.size(), timer.elapsed());

    return coloring;
}

//! Traits specifying if a given discretization tag supports coloring
template<class DiscretizationMethod>
struct SupportsColoring : public std::false_type {};

template<> struct SupportsColoring<DiscretizationMethods::CCTpfa> : public std::true_type {};
template<> struct SupportsColoring<DiscretizationMethods::CCMpfa> : public std::true_type {};
template<> struct SupportsColoring<DiscretizationMethods::Box> : public std::true_type {};
template<> struct SupportsColoring<DiscretizationMethods::FCStaggered> : public std::true_type {};
template<> struct SupportsColoring<DiscretizationMethods::FCDiamond> : public std::true_type {};
template<> struct SupportsColoring<DiscretizationMethods::PQ1Bubble> : public std::true_type {};

} // end namespace Dumux

#endif
