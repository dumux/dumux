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
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Geometry
 * \brief Algorithms that finds which geometric entites intersect
 */
#ifndef DUMUX_INTERSECTING_ENTITIES_HH
#define DUMUX_INTERSECTING_ENTITIES_HH

#include <cmath>
#include <type_traits>

#include <dune/common/fvector.hh>

#include <dumux/common/math.hh>
#include <dumux/common/geometry/boundingboxtree.hh>
#include <dumux/common/geometry/intersectspointgeometry.hh>
#include <dumux/common/geometry/geometryintersection.hh>
#include <dumux/common/geometry/boundingboxtreeintersection.hh>
#include <dumux/common/geometry/triangulation.hh>

namespace Dumux {

/*!
 * \ingroup Geometry
 * \brief Compute all intersections between entities and a point
 */
template<class EntitySet, class ctype, int dimworld>
inline std::vector<std::size_t>
intersectingEntities(const Dune::FieldVector<ctype, dimworld>& point,
                     const BoundingBoxTree<EntitySet>& tree,
                     bool isCartesianGrid = false)
{
    // Call the recursive find function to find candidates
    std::vector<std::size_t> entities;
    intersectingEntities(point, tree, tree.numBoundingBoxes() - 1, entities, isCartesianGrid);
    return entities;
}

/*!
 * \ingroup Geometry
 * \brief Compute intersections with point for all nodes of the bounding box tree recursively
 */
template<class EntitySet, class ctype, int dimworld>
void intersectingEntities(const Dune::FieldVector<ctype, dimworld>& point,
                          const BoundingBoxTree<EntitySet>& tree,
                          std::size_t node,
                          std::vector<std::size_t>& entities,
                          bool isCartesianGrid = false)
{
    // Get the bounding box for the current node
    const auto& bBox = tree.getBoundingBoxNode(node);

    // if the point is not in the bounding box we can stop
    if (!intersectsPointBoundingBox(point, tree.getBoundingBoxCoordinates(node))) return;

    // now we know it's inside
    // if the box is a leaf do the primitive test.
    else if (tree.isLeaf(bBox, node))
    {
        const std::size_t entityIdx = bBox.child1;
        // for structured cube grids skip the primitive test
        if (isCartesianGrid)
            entities.push_back(entityIdx);
        else
        {
            const auto geometry = tree.entitySet().entity(entityIdx).geometry();
            // if the primitive is positive it intersects the actual geometry, add the entity to the list
            if (intersectsPointGeometry(point, geometry))
                entities.push_back(entityIdx);
        }
    }

    // No leaf. Check both children nodes.
    else
    {
        intersectingEntities(point, tree, bBox.child0, entities, isCartesianGrid);
        intersectingEntities(point, tree, bBox.child1, entities, isCartesianGrid);
    }
}

/*!
 * \ingroup Geometry
 * \brief Compute all intersections between two bounding box trees
 */
template<class EntitySet0, class EntitySet1>
inline std::vector<BoundingBoxTreeIntersection<EntitySet0, EntitySet1>>
intersectingEntities(const BoundingBoxTree<EntitySet0>& treeA,
                     const BoundingBoxTree<EntitySet1>& treeB)
{
    // check if the world dimensions match
    static_assert(int(EntitySet0::dimensionworld) == int(EntitySet1::dimensionworld),
        "Can only intersect bounding box trees of same world dimension");

    // Create data structure for return type
    std::vector<BoundingBoxTreeIntersection<EntitySet0, EntitySet1>> intersections;

    // Call the recursive find function to find candidates
    intersectingEntities(treeA, treeB,
                         treeA.numBoundingBoxes() - 1,
                         treeB.numBoundingBoxes() - 1,
                         intersections);

    return intersections;
}

/*!
 * \ingroup Geometry
 * \brief Compute all intersections between two all bounding box tree nodes recursively
 */
template<class EntitySet0, class EntitySet1>
void intersectingEntities(const BoundingBoxTree<EntitySet0>& treeA,
                          const BoundingBoxTree<EntitySet1>& treeB,
                          std::size_t nodeA, std::size_t nodeB,
                          std::vector<BoundingBoxTreeIntersection<EntitySet0, EntitySet1>>& intersections)
{
    // Get the bounding box for the current node
    const auto& bBoxA = treeA.getBoundingBoxNode(nodeA);
    const auto& bBoxB = treeB.getBoundingBoxNode(nodeB);

    // if the two bounding boxes of the current nodes don't intersect we can stop searching
    static constexpr int dimworld = EntitySet0::dimensionworld;
    if (!intersectsBoundingBoxBoundingBox<dimworld>(treeA.getBoundingBoxCoordinates(nodeA),
                                                    treeB.getBoundingBoxCoordinates(nodeB)))
        return;

    // Check if we have a leaf in treeA or treeB
    const bool isLeafA = treeA.isLeaf(bBoxA, nodeA);
    const bool isLeafB = treeB.isLeaf(bBoxB, nodeB);

    // If both boxes are leaves do the primitive test
    if (isLeafA && isLeafB)
    {
        const auto eIdxA = bBoxA.child1;
        const auto eIdxB = bBoxB.child1;

        const auto geometryA = treeA.entitySet().entity(eIdxA).geometry();
        const auto geometryB = treeB.entitySet().entity(eIdxB).geometry();

        using GeometryA = std::decay_t<decltype(geometryA)>;
        using GeometryB = std::decay_t<decltype(geometryB)>;
        using Policy = IntersectionPolicy::DefaultPolicy<GeometryA, GeometryB>;
        using IntersectionAlgorithm = GeometryIntersection<GeometryA, GeometryB, Policy>;
        using Intersection = typename IntersectionAlgorithm::Intersection;
        Intersection intersection;

        if (IntersectionAlgorithm::intersection(geometryA, geometryB, intersection))
        {
            static constexpr int dimIntersection = Policy::dimIntersection;

            if (dimIntersection >= 2)
            {
                const auto triangulation = triangulate<dimIntersection, dimworld>(intersection);
                for (unsigned int i = 0; i < triangulation.size(); ++i)
                    intersections.emplace_back(eIdxA, eIdxB, std::move(triangulation[i]));
            }
            else
                intersections.emplace_back(eIdxA, eIdxB, intersection);
        }
    }

    // if we reached the leaf in treeA, just continue in treeB
    else if (isLeafA)
    {
        intersectingEntities(treeA, treeB, nodeA, bBoxB.child0, intersections);
        intersectingEntities(treeA, treeB, nodeA, bBoxB.child1, intersections);
    }

    // if we reached the leaf in treeB, just continue in treeA
    else if (isLeafB)
    {
        intersectingEntities(treeA, treeB, bBoxA.child0, nodeB, intersections);
        intersectingEntities(treeA, treeB, bBoxA.child1, nodeB, intersections);
    }

    // we know now that both trees didn't reach the leaf yet so
    // we continue with the larger tree first (bigger node number)
    else if (nodeA > nodeB)
    {
        intersectingEntities(treeA, treeB, bBoxA.child0, nodeB, intersections);
        intersectingEntities(treeA, treeB, bBoxA.child1, nodeB, intersections);
    }
    else
    {
        intersectingEntities(treeA, treeB, nodeA, bBoxB.child0, intersections);
        intersectingEntities(treeA, treeB, nodeA, bBoxB.child1, intersections);
    }
}

} // end namespace Dumux

#endif
