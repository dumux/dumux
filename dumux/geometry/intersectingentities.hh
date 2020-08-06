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
#ifndef DUMUX_GEOMETRY_INTERSECTING_ENTITIES_HH
#define DUMUX_GEOMETRY_INTERSECTING_ENTITIES_HH

#include <cmath>
#include <type_traits>
#include <vector>
#include <algorithm>

#include <dune/common/fvector.hh>

#include <dumux/common/math.hh>
#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/geometry/intersectspointgeometry.hh>
#include <dumux/geometry/geometryintersection.hh>
#include <dumux/geometry/triangulation.hh>

namespace Dumux {

/*!
 * \ingroup Geometry
 * \brief An intersection object resulting from the intersection of two primitives in an entity set
 * \note this is used as return type for some of the intersectingEntities overloads below
 */
template<int dimworld, class CoordTypeA, class CoordTypeB = CoordTypeA>
class IntersectionInfo
{
public:
    using ctype = typename Dune::PromotionTraits<CoordTypeA, CoordTypeB>::PromotedType;
    static constexpr int dimensionworld = dimworld;
    using GlobalPosition = Dune::FieldVector<ctype, dimworld>;

    template<class Corners>
    explicit IntersectionInfo(std::size_t a, std::size_t b, Corners&& c)
    : a_(a)
    , b_(b)
    , corners_(c.begin(), c.end())
    {}

    //! Get the index of the intersecting entity belonging to this grid
    std::size_t first() const
    { return a_; }

    //! Get the index of the intersecting entity belonging to the other grid
    std::size_t second() const
    { return b_; }

    //! Get the corners of the intersection geometry
    std::vector<GlobalPosition> corners() const
    { return corners_; }

    /*!
     * \brief Check if the corners of this intersection match with the given corners
     * \note This is useful to check if the intersection geometry of two intersections coincide.
     */
    bool cornersMatch(const std::vector<GlobalPosition>& otherCorners) const
    {
        if (otherCorners.size() != corners_.size())
            return false;

        const auto eps = 1.5e-7*(corners_[1] - corners_[0]).two_norm();
        for (int i = 0; i < corners_.size(); ++i)
            if ((corners_[i] - otherCorners[i]).two_norm() > eps)
                return false;

        return true;
    }

private:
    std::size_t a_, b_; //!< Indices of the intersection elements
    std::vector<GlobalPosition> corners_; //!< the corner points of the intersection geometry
};

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
 * \brief Compute all intersections between a geometry and a bounding box tree
 */
template<class Geometry, class EntitySet>
inline std::vector<IntersectionInfo<Geometry::coorddimension, typename Geometry::ctype, typename EntitySet::ctype>>
intersectingEntities(const Geometry& geometry,
                     const BoundingBoxTree<EntitySet>& tree)
{
    // check if the world dimensions match
    static_assert(int(Geometry::coorddimension) == int(EntitySet::dimensionworld),
        "Can only intersect geometry and bounding box tree of same world dimension");

    // Create data structure for return type
    std::vector<IntersectionInfo<Geometry::coorddimension, typename Geometry::ctype, typename EntitySet::ctype>> intersections;
    using ctype = typename IntersectionInfo<Geometry::coorddimension, typename Geometry::ctype, typename EntitySet::ctype>::ctype;
    static constexpr int dimworld = Geometry::coorddimension;

    // compute the bounding box of the given geometry
    std::array<ctype, 2*Geometry::coorddimension> bBox;
    ctype* xMin = bBox.data(); ctype* xMax = xMin + Geometry::coorddimension;

    // Get coordinates of first vertex
    auto corner = geometry.corner(0);
    for (std::size_t dimIdx = 0; dimIdx < dimworld; ++dimIdx)
        xMin[dimIdx] = xMax[dimIdx] = corner[dimIdx];

    // Compute the min and max over the remaining vertices
    for (std::size_t cornerIdx = 1; cornerIdx < geometry.corners(); ++cornerIdx)
    {
        corner = geometry.corner(cornerIdx);
        for (std::size_t dimIdx = 0; dimIdx < dimworld; ++dimIdx)
        {
            using std::max;
            using std::min;
            xMin[dimIdx] = min(xMin[dimIdx], corner[dimIdx]);
            xMax[dimIdx] = max(xMax[dimIdx], corner[dimIdx]);
        }
    }

    // Call the recursive find function to find candidates
    intersectingEntities(geometry, tree,
                         bBox, tree.numBoundingBoxes() - 1,
                         intersections);

    return intersections;
}

/*!
 * \ingroup Geometry
 * \brief Compute intersections with point for all nodes of the bounding box tree recursively
 */
template<class Geometry, class EntitySet>
void intersectingEntities(const Geometry& geometry,
                          const BoundingBoxTree<EntitySet>& tree,
                          const std::array<typename Geometry::ctype, 2*Geometry::coorddimension>& bBox,
                          std::size_t nodeIdx,
                          std::vector<IntersectionInfo<Geometry::coorddimension,
                                                       typename Geometry::ctype,
                                                       typename EntitySet::ctype>>& intersections)
{
    // if the two bounding boxes don't intersect we can stop searching
    static constexpr int dimworld = Geometry::coorddimension;
    if (!intersectsBoundingBoxBoundingBox<dimworld>(bBox.data(), tree.getBoundingBoxCoordinates(nodeIdx)))
        return;

    // get node info for current bounding box node
    const auto& bBoxNode = tree.getBoundingBoxNode(nodeIdx);

    // if the box is a leaf do the primitive test.
    if (tree.isLeaf(bBoxNode, nodeIdx))
    {
        // eIdxA is always 0 since we intersect with exactly one geometry
        const auto eIdxA = 0;
        const auto eIdxB = bBoxNode.child1;

        const auto geometryTree = tree.entitySet().entity(eIdxB).geometry();
        using GeometryTree = std::decay_t<decltype(geometryTree)>;
        using Policy = IntersectionPolicy::DefaultPolicy<Geometry, GeometryTree>;
        using IntersectionAlgorithm = GeometryIntersection<Geometry, GeometryTree, Policy>;
        using Intersection = typename IntersectionAlgorithm::Intersection;
        Intersection intersection;

        if (IntersectionAlgorithm::intersection(geometry, geometryTree, intersection))
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

    // No leaf. Check both children nodes.
    else
    {
        intersectingEntities(geometry, tree, bBox, bBoxNode.child0, intersections);
        intersectingEntities(geometry, tree, bBox, bBoxNode.child1, intersections);
    }
}

/*!
 * \ingroup Geometry
 * \brief Compute all intersections between two bounding box trees
 */
template<class EntitySet0, class EntitySet1>
inline std::vector<IntersectionInfo<EntitySet0::dimensionworld, typename EntitySet0::ctype, typename EntitySet1::ctype>>
intersectingEntities(const BoundingBoxTree<EntitySet0>& treeA,
                     const BoundingBoxTree<EntitySet1>& treeB)
{
    // check if the world dimensions match
    static_assert(int(EntitySet0::dimensionworld) == int(EntitySet1::dimensionworld),
        "Can only intersect bounding box trees of same world dimension");

    // Create data structure for return type
    std::vector<IntersectionInfo<EntitySet0::dimensionworld, typename EntitySet0::ctype, typename EntitySet1::ctype>> intersections;

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
                          std::vector<IntersectionInfo<EntitySet0::dimensionworld,
                                                       typename EntitySet0::ctype,
                                                       typename EntitySet1::ctype>>& intersections)
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

/*!
 * \ingroup Geometry
 * \brief Compute the index of the intersecting element of a Cartesian grid with a point
 * The grid is given by the lower left corner (min), the upper right corner (max)
 * and the number of cells in each direction (cells).
 * \note If there are several options the lowest matching cell index will be returned
 * \note Returns the index i + I*j + I*J*k for the intersecting element i,j,k of a grid with I,J,K cells
 */
template<class ctype, int dimworld>
inline std::size_t intersectingEntityCartesianGrid(const Dune::FieldVector<ctype, dimworld>& point,
                                                   const Dune::FieldVector<ctype, dimworld>& min,
                                                   const Dune::FieldVector<ctype, dimworld>& max,
                                                   const std::array<int, std::size_t(dimworld)>& cells)
{
    std::size_t index = 0;
    for (int i = 0; i < dimworld; ++i)
    {
        using std::clamp; using std::floor;
        ctype dimOffset = clamp<ctype>(floor((point[i]-min[i])*cells[i]/(max[i]-min[i])), 0.0, cells[i]-1);
        for (int j = 0; j < i; ++j)
            dimOffset *= cells[j];
        index += static_cast<std::size_t>(dimOffset);
    }
    return index;
}

} // end namespace Dumux

#endif
