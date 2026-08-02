// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Nonlinear
 * \brief Building subdomain decompositions from a grid geometry
 */
#ifndef DUMUX_NONLINEAR_PRECONDITIONING_PARTITIONER_HH
#define DUMUX_NONLINEAR_PRECONDITIONING_PARTITIONER_HH

#include <algorithm>
#include <array>
#include <cstddef>
#include <memory>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dumux/nonlinear/preconditioning/partition.hh>

namespace Dumux::NonlinearPreconditioning {

/*!
 * \ingroup Nonlinear
 * \brief The influence map of a cell-centred scheme: which rows a visited element contributes to
 *
 * This is the connectivity map the assembler itself uses, so subdomain boundaries built from it always
 * align with true stencil boundaries.
 */
template<class GridGeometry>
std::vector<std::vector<Partition::Index>> buildInfluenceMap(const GridGeometry& gridGeometry)
{
    const auto numElements = gridGeometry.gridView().size(0);
    std::vector<std::vector<Partition::Index>> influence(numElements);
    for (Partition::Index i = 0; i < static_cast<Partition::Index>(numElements); ++i)
        for (const auto& dataJ : gridGeometry.connectivityMap()[i])
            influence[i].push_back(dataJ.globalJ);
    return influence;
}

/*!
 * \ingroup Nonlinear
 * \brief Decompose the bounding box into an axis-aligned block grid of subdomains
 * \param gridGeometry the grid geometry
 * \param influence the influence map, see buildInfluenceMap()
 * \param blocksPerDirection number of subdomains along each coordinate direction
 * \param overlapWidth number of connectivity rings the cores are grown by
 *
 * Cheap, deterministic and independent of any graph partitioning library, which makes it the sensible
 * default on structured grids. Cells are binned by the position of their centre, so a block may come out
 * empty on an unstructured or locally refined grid; such blocks are dropped rather than left to produce
 * an empty subdomain.
 */
template<class GridGeometry>
std::shared_ptr<Partition> makeCartesianPartition(
    const GridGeometry& gridGeometry,
    const std::vector<std::vector<Partition::Index>>& influence,
    const std::array<std::size_t, GridGeometry::GridView::dimensionworld>& blocksPerDirection,
    std::size_t overlapWidth = 0)
{
    static constexpr int dimWorld = GridGeometry::GridView::dimensionworld;

    std::size_t numBlocks = 1;
    for (const auto n : blocksPerDirection)
    {
        if (n == 0)
            DUNE_THROW(Dune::InvalidStateException, "Number of blocks per direction must be positive");
        numBlocks *= n;
    }

    const auto& bBoxMin = gridGeometry.bBoxMin();
    const auto& bBoxMax = gridGeometry.bBoxMax();

    std::vector<SubdomainIndex> raw(gridGeometry.gridView().size(0), noSubdomain);
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        const auto eIdx = gridGeometry.elementMapper().index(element);
        const auto center = element.geometry().center();

        std::size_t block = 0;
        for (int dir = 0; dir < dimWorld; ++dir)
        {
            const auto extent = bBoxMax[dir] - bBoxMin[dir];
            const auto relative = extent > 0.0 ? (center[dir] - bBoxMin[dir])/extent : 0.0;
            auto k = static_cast<std::size_t>(relative*double(blocksPerDirection[dir]));
            k = std::min(k, blocksPerDirection[dir] - 1);
            block = block*blocksPerDirection[dir] + k;
        }
        raw[eIdx] = block;
    }

    std::vector<SubdomainIndex> compressed(numBlocks, noSubdomain);
    SubdomainIndex next = 0;
    for (const auto block : raw)
        if (compressed[block] == noSubdomain)
            compressed[block] = next++;

    std::vector<SubdomainIndex> cores(raw.size());
    for (std::size_t i = 0; i < raw.size(); ++i)
        cores[i] = compressed[raw[i]];

    return std::make_shared<Partition>(influence, std::move(cores), overlapWidth);
}

/*!
 * \ingroup Nonlinear
 * \brief Grow subdomains from seeds by breadth-first search over the connectivity graph
 * \note Dependency-free alternative to a graph partitioner, for grids where binning by position is a
 *       poor fit. Seeds are spread evenly over the element index range.
 */
inline std::shared_ptr<Partition> makeGreedyPartition(
    const std::vector<std::vector<Partition::Index>>& influence,
    std::size_t numSubdomains,
    std::size_t overlapWidth = 0)
{
    using Index = Partition::Index;
    const auto numElements = influence.size();
    if (numSubdomains == 0 || numSubdomains > numElements)
        DUNE_THROW(Dune::InvalidStateException, "Invalid number of subdomains " << numSubdomains);

    std::vector<std::vector<Index>> neighbours(numElements);
    for (Index i = 0; i < numElements; ++i)
        for (const auto j : influence[i])
        {
            neighbours[i].push_back(j);
            neighbours[j].push_back(i);
        }

    std::vector<SubdomainIndex> cores(numElements, noSubdomain);
    std::vector<std::vector<Index>> fronts(numSubdomains);
    for (SubdomainIndex s = 0; s < numSubdomains; ++s)
    {
        const auto seed = static_cast<Index>(s*numElements/numSubdomains);
        cores[seed] = s;
        fronts[s].push_back(seed);
    }

    // grown one ring at a time in turn, so the subdomains stay comparable in size
    bool grew = true;
    while (grew)
    {
        grew = false;
        for (SubdomainIndex s = 0; s < numSubdomains; ++s)
        {
            std::vector<Index> next;
            for (const auto element : fronts[s])
                for (const auto j : neighbours[element])
                    if (cores[j] == noSubdomain)
                    {
                        cores[j] = s;
                        next.push_back(j);
                    }
            if (!next.empty())
                grew = true;
            fronts[s] = std::move(next);
        }
    }

    // a disconnected graph can leave elements unclaimed; give them to their lowest-indexed neighbour
    for (Index i = 0; i < numElements; ++i)
        if (cores[i] == noSubdomain)
            cores[i] = 0;

    return std::make_shared<Partition>(influence, std::move(cores), overlapWidth);
}

} // end namespace Dumux::NonlinearPreconditioning

#endif
