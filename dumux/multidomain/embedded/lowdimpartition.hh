// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedCoupling
 * \brief Partition of an embedded lower-dimensional grid that follows a distributed bulk grid
 */
#ifndef DUMUX_MULTIDOMAIN_EMBEDDED_LOWDIM_PARTITION_HH
#define DUMUX_MULTIDOMAIN_EMBEDDED_LOWDIM_PARTITION_HH

#include <cstddef>
#include <limits>
#include <memory>
#include <type_traits>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/grid/common/partitionset.hh>

#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/geometry/distance.hh>
#include <dumux/geometry/geometricentityset.hh>
#include <dumux/geometry/intersectingentities.hh>

namespace Dumux {

/*!
 * \ingroup EmbeddedCoupling
 * \brief The rank each element of an embedded lower-dimensional grid belongs to, given a
 *        distributed bulk grid
 *
 * An element belongs to the rank whose interior bulk elements contain its center. If no bulk element
 * contains the center, e.g. because the lower-dimensional grid lies in a hole of the bulk grid or
 * outside of it, the element belongs to the rank with the closest bulk domain boundary. Ties go to
 * the lowest rank.
 *
 * Collective over the bulk grid's communicator. Every rank has to hold the whole lower-dimensional
 * grid and gets the same result.
 *
 * \param bulkGridView the bulk grid view, distributed
 * \param lowDimGridView the lower-dimensional grid view, the same on every rank
 * \return the rank per lower-dimensional element, indexed by the element index of lowDimGridView
 */
template<class BulkGridView, class LowDimGridView>
std::vector<int> lowDimPartition(const BulkGridView& bulkGridView, const LowDimGridView& lowDimGridView)
{
    using Scalar = typename BulkGridView::ctype;
    static constexpr int dim = BulkGridView::dimension;
    static constexpr int dimWorld = BulkGridView::dimensionworld;

    // geometries of some grids refer to storage of the grid that is reused, so keep copies
    const auto copy = [](const auto& geometry)
    {
        using Geometry = std::decay_t<decltype(geometry)>;
        std::vector<typename Geometry::GlobalCoordinate> corners(geometry.corners());
        for (int i = 0; i < geometry.corners(); ++i)
            corners[i] = geometry.corner(i);
        return Dune::MultiLinearGeometry<Scalar, Geometry::mydimension, dimWorld>(geometry.type(), corners);
    };

    using ElementGeometry = Dune::MultiLinearGeometry<Scalar, dim, dimWorld>;
    using FacetGeometry = Dune::MultiLinearGeometry<Scalar, dim-1, dimWorld>;
    std::vector<ElementGeometry> elementGeometries;
    std::vector<FacetGeometry> boundaryFacetGeometries;
    for (const auto& element : elements(bulkGridView, Dune::Partitions::interior))
    {
        elementGeometries.push_back(copy(element.geometry()));
        for (const auto& intersection : intersections(bulkGridView, element))
            if (intersection.boundary())
                boundaryFacetGeometries.push_back(copy(intersection.geometry()));
    }

    using ElementSet = GeometriesEntitySet<ElementGeometry>;
    using FacetSet = GeometriesEntitySet<FacetGeometry>;
    std::unique_ptr<BoundingBoxTree<ElementSet>> elementTree;
    std::unique_ptr<BoundingBoxTree<FacetSet>> facetTree;
    if (!elementGeometries.empty())
        elementTree = std::make_unique<BoundingBoxTree<ElementSet>>(std::make_shared<ElementSet>(std::move(elementGeometries)));
    if (!boundaryFacetGeometries.empty())
        facetTree = std::make_unique<BoundingBoxTree<FacetSet>>(std::make_shared<FacetSet>(std::move(boundaryFacetGeometries)));

    // the lowest rank whose interior bulk elements contain the center
    const auto& comm = bulkGridView.comm();
    const auto& indexSet = lowDimGridView.indexSet();
    std::vector<typename LowDimGridView::template Codim<0>::Geometry::GlobalCoordinate> centers(lowDimGridView.size(0));
    std::vector<int> rank(centers.size(), comm.size());
    for (const auto& element : elements(lowDimGridView))
    {
        const auto idx = indexSet.index(element);
        centers[idx] = element.geometry().center();
        if (elementTree && !intersectingEntities(centers[idx], *elementTree).empty())
            rank[idx] = comm.rank();
    }
    if (!rank.empty())
        comm.min(rank.data(), rank.size());

    // the elements whose center no rank contains go to the rank with the closest bulk domain boundary
    std::vector<std::size_t> outside;
    for (std::size_t i = 0; i < rank.size(); ++i)
        if (rank[i] == comm.size())
            outside.push_back(i);

    if (!outside.empty())
    {
        constexpr Scalar unreachable = std::numeric_limits<Scalar>::max();
        std::vector<Scalar> squaredDistance(outside.size(), unreachable);
        if (facetTree)
            for (std::size_t k = 0; k < outside.size(); ++k)
                squaredDistance[k] = closestEntity(centers[outside[k]], *facetTree).first;

        auto minSquaredDistance = squaredDistance;
        comm.min(minSquaredDistance.data(), minSquaredDistance.size());

        std::vector<int> outsideRank(outside.size(), comm.size());
        for (std::size_t k = 0; k < outside.size(); ++k)
            if (squaredDistance[k] == minSquaredDistance[k] && squaredDistance[k] < unreachable)
                outsideRank[k] = comm.rank();
        comm.min(outsideRank.data(), outsideRank.size());

        for (std::size_t k = 0; k < outside.size(); ++k)
            rank[outside[k]] = outsideRank[k];
    }

    for (std::size_t i = 0; i < rank.size(); ++i)
        if (rank[i] == comm.size())
            DUNE_THROW(Dune::InvalidStateException, "No rank holds a bulk element or boundary facet for lower-dimensional element " << i);

    return rank;
}

} // end namespace Dumux

#endif
