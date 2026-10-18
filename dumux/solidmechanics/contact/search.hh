// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Contact
 * \brief Search for the facets near the nodes of a contact surface
 */
#ifndef DUMUX_SOLIDMECHANICS_CONTACT_SEARCH_HH
#define DUMUX_SOLIDMECHANICS_CONTACT_SEARCH_HH

#include <algorithm>
#include <memory>
#include <vector>

#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/type.hh>

#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/geometry/geometricentityset.hh>
#include <dumux/solidmechanics/contact/closestpoint.hh>

namespace Dumux::Contact {

namespace Detail {

template<class Tree, class Point, class Scalar>
void collectBoxesWithinDistance(const Tree& tree, std::size_t nodeIdx, const Point& point,
                                Scalar squaredRadius, std::vector<std::size_t>& entities)
{
    if (squaredDistancePointBoundingBox(point, tree.getBoundingBoxCoordinates(nodeIdx)) > squaredRadius)
        return;

    const auto& node = tree.getBoundingBoxNode(nodeIdx);
    if (tree.isLeaf(node, nodeIdx))
        entities.push_back(node.child1);
    else
    {
        collectBoxesWithinDistance(tree, node.child0, point, squaredRadius, entities);
        collectBoxesWithinDistance(tree, node.child1, point, squaredRadius, entities);
    }
}

} // end namespace Detail

/*!
 * \ingroup Contact
 * \brief For every surface node, the facets within a distance and accepted by a predicate
 *
 * A bounding box tree over the facets in the given configuration selects the facets whose box is
 * within the radius; the exact distance to the facet decides.
 *
 * \param surface the contact surface
 * \param x the positions of the surface nodes
 * \param radius the search radius
 * \param accept returns whether a facet is a candidate for a node, called as accept(node, facet)
 * \return the sorted candidate facets of every surface node
 */
template<class Surface, class Scalar, class Accept>
std::vector<std::vector<std::size_t>>
findCandidateFacets(const Surface& surface,
                    const std::vector<typename Surface::GlobalPosition>& x,
                    Scalar radius, const Accept& accept)
{
    static constexpr int dim = Surface::dimension;
    using Geometry = Dune::MultiLinearGeometry<Scalar, dim-1, dim>;
    std::vector<Geometry> geometries;
    geometries.reserve(surface.numFacets());
    for (std::size_t f = 0; f < surface.numFacets(); ++f)
    {
        const auto c = surface.corners(f, x);
        geometries.emplace_back(Dune::GeometryTypes::simplex(dim-1), std::vector(c.begin(), c.end()));
    }
    using EntitySet = GeometriesEntitySet<Geometry>;
    const BoundingBoxTree<EntitySet> tree(std::make_shared<EntitySet>(std::move(geometries)));

    std::vector<std::vector<std::size_t>> candidates(surface.numNodes());
    std::vector<std::size_t> boxes;
    const Scalar squaredRadius = radius*radius;
    for (std::size_t s = 0; s < surface.numNodes(); ++s)
    {
        boxes.clear();
        Detail::collectBoxesWithinDistance(tree, tree.numBoundingBoxes() - 1, x[s], squaredRadius, boxes);
        for (const auto f : boxes)
            if (closestPointOnSimplex(x[s], surface.corners(f, x)).squaredDistance <= squaredRadius && accept(s, f))
                candidates[s].push_back(f);
        std::sort(candidates[s].begin(), candidates[s].end());
    }
    return candidates;
}

/*!
 * \ingroup Contact
 * \brief The same as findCandidateFacets by testing all pairs of nodes and facets
 */
template<class Surface, class Scalar, class Accept>
std::vector<std::vector<std::size_t>>
findCandidateFacetsBruteForce(const Surface& surface,
                              const std::vector<typename Surface::GlobalPosition>& x,
                              Scalar radius, const Accept& accept)
{
    std::vector<std::vector<std::size_t>> candidates(surface.numNodes());
    for (std::size_t s = 0; s < surface.numNodes(); ++s)
        for (std::size_t f = 0; f < surface.numFacets(); ++f)
            if (closestPointOnSimplex(x[s], surface.corners(f, x)).squaredDistance <= radius*radius && accept(s, f))
                candidates[s].push_back(f);
    return candidates;
}

} // end namespace Dumux::Contact

#endif
