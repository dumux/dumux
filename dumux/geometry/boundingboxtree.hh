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
 * \brief An axis-aligned bounding box volume hierarchy for dune grids
 *
 * Dumux implementation of an AABB tree
 * Inspired by the AABB tree implementation in DOLFIN by Anders Logg which has the
 * following license info: DOLFIN is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 */
#ifndef DUMUX_GEOMETRY_BOUNDINGBOXTREE_HH
#define DUMUX_GEOMETRY_BOUNDINGBOXTREE_HH

#include <vector>
#include <array>
#include <algorithm>
#include <memory>
#include <numeric>
#include <type_traits>
#include <iostream>

#include <dune/common/promotiontraits.hh>
#include <dune/common/timer.hh>
#include <dune/common/fvector.hh>

namespace Dumux {

/*!
 * \ingroup Geometry
 * \brief An axis-aligned bounding box volume tree implementation
 *
 * The class constructs a hierarchical structure of bounding box volumes around
 * grid entities. This class can be used to efficiently compute intersections
 * between a grid and other geometrical object. It only implements the intersection
 * of two of such bounding box trees, so that two independent grids can be intersected.
 * \tparam GeometricEntitySet has the following requirements
 *         * export dimensionworld, ctype
 *         * a size() member function returning the number of entities
 *         * begin() and end() member function returning at least forward iterators to entities
 *         * an index() method returning a consecutive index given an entity
 *         * an entity() method returning an entity given the consecutive index
 *         * entities have the following requirements:
 *             * a member function geometry() returning a geometry with the member functions
 *                 * corner() and corners() returning global coordinates and number of corners
 */
template <class GeometricEntitySet>
class BoundingBoxTree
{
    enum { dimworld = GeometricEntitySet::dimensionworld };
    using ctype = typename GeometricEntitySet::ctype;

    /*!
     * \brief Bounding box node data structure
     * leaf nodes are indicated by setting child0 to
     * the node itself and child1 to the index of the entity in the bounding box.
     */
    struct BoundingBoxNode
    {
        std::size_t child0;
        std::size_t child1;
    };

public:
    //! the type of entity set this tree was built with
    using EntitySet = GeometricEntitySet;

    //! Default Constructor
    BoundingBoxTree() = default;

    //! Constructor with gridView
    BoundingBoxTree(std::shared_ptr<const GeometricEntitySet> set)
    { build(set); }

    //! Build up bounding box tree for a grid with leafGridView
    void build(std::shared_ptr<const GeometricEntitySet> set)
    {
        // set the pointer to the entity set
        entitySet_ = set;

        // clear all internal data
        boundingBoxNodes_.clear();
        boundingBoxCoordinates_.clear();

        // start the timer
        Dune::Timer timer;

        // Create bounding boxes for all elements
        const auto numLeaves = set->size();

        // reserve enough space for the nodes and the coordinates
        const auto numNodes = 2*numLeaves - 1;
        boundingBoxNodes_.reserve(numNodes);
        boundingBoxCoordinates_.reserve(numNodes*2*dimworld);

        // create a vector for leaf boxes (min and max for all dims)
        std::vector<ctype> leafBoxes(2*dimworld*numLeaves);

        for (const auto& geometricEntity : *set)
            computeEntityBoundingBox_(leafBoxes.data() + 2*dimworld*set->index(geometricEntity), geometricEntity);

        // create the leaf partition, the set of available indices (to be sorted)
        std::vector<std::size_t> leafPartition(numLeaves);
        std::iota(leafPartition.begin(), leafPartition.end(), 0);

        // Recursively build the bounding box tree
        build_(leafBoxes, leafPartition.begin(), leafPartition.end());

        // We are done, log output
        std::cout << "Computed bounding box tree with " << numBoundingBoxes()
                  << " nodes for " << numLeaves << " grid entites in "
                  << timer.stop() << " seconds." << std::endl;
    }

    //! the entity set this tree was built with
    const EntitySet& entitySet() const
    { return *entitySet_; }

    /////////////////////////////////////////////////////
    //! Interface to be used by other bounding box trees
    /////////////////////////////////////////////////////

    //! Get an existing bounding box for a given node
    const BoundingBoxNode& getBoundingBoxNode(std::size_t nodeIdx) const
    { return boundingBoxNodes_[nodeIdx]; }

    //! Get an existing bounding box for a given node
    const ctype* getBoundingBoxCoordinates(std::size_t nodeIdx) const
    { return boundingBoxCoordinates_.data() + 2*dimworld*nodeIdx; }

    //! Get the number of bounding boxes currently in the tree
    std::size_t numBoundingBoxes() const
    { return boundingBoxNodes_.size(); }

    //! Check whether a bounding box node is a leaf node
    //! Leaf nodes have itself as child0
    bool isLeaf(const BoundingBoxNode& node, std::size_t nodeIdx) const
    { return node.child0 == nodeIdx; }

private:

    //! vector of bounding box nodes
    std::vector<BoundingBoxNode> boundingBoxNodes_;

    //! vector of bounding box coordinates
    std::vector<ctype> boundingBoxCoordinates_;

    //! a pointer to the entity set
    std::shared_ptr<const EntitySet> entitySet_;

    //! Compute the bounding box of a grid entity
    template <class Entity>
    void computeEntityBoundingBox_(ctype* b, const Entity& entity) const
    {
        // get the bounding box coordinates
        ctype* xMin = b;
        ctype* xMax = b + dimworld;

        // get mesh entity data
        auto geometry = entity.geometry();

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
    }

    //! Build bounding box tree for all entities recursively
    std::size_t build_(const std::vector<ctype>& leafBoxes,
                       const std::vector<std::size_t>::iterator& begin,
                       const std::vector<std::size_t>::iterator& end)
    {
        assert(begin < end);

        // If we reached the end of the recursion, i.e. only a leaf box is left
        if (end - begin == 1)
        {
            // Get the bounding box coordinates for the leaf
            const std::size_t leafNodeIdx = *begin;
            const auto beginCoords = leafBoxes.begin() + 2*dimworld*leafNodeIdx;
            const auto endCoords = beginCoords + 2*dimworld;

            // Store the data in the bounding box
            // leaf nodes are indicated by setting child0 to
            // the node itself and child1 to the index of the entity in the bounding box.
            return addBoundingBox_(BoundingBoxNode{numBoundingBoxes(), leafNodeIdx}, beginCoords, endCoords);
        }

        // Compute the bounding box of all bounding boxes in the range [begin, end]
        const auto bCoords = computeBBoxOfBBoxes_(leafBoxes, begin, end);

        // sort bounding boxes along the longest axis
        const auto axis = computeLongestAxis_(bCoords);

        // nth_element sorts the range to make sure that middle points to the coordinate median in axis direction
        // this is the most expensive part of the algorithm
        auto middle = begin + (end - begin)/2;
        std::nth_element(begin, middle, end, [&leafBoxes, axis](std::size_t i, std::size_t j)
                         {
                             const ctype* bi = leafBoxes.data() + 2*dimworld*i;
                             const ctype* bj = leafBoxes.data() + 2*dimworld*j;
                             return bi[axis] + bi[axis + dimworld] < bj[axis] + bj[axis + dimworld];
                         });

        // split the bounding boxes into two at the middle iterator and call build recursively, each
        // call resulting in a new node of this bounding box, i.e. the root will be added at the end of the process.
        return addBoundingBox_(BoundingBoxNode{build_(leafBoxes, begin, middle), build_(leafBoxes, middle, end)},
                               bCoords.begin(), bCoords.end());
    }

    //! Add a new bounding box to the tree
    template <class Iterator>
    std::size_t addBoundingBox_(BoundingBoxNode&& node,
                                const Iterator& coordBegin,
                                const Iterator& coordEnd)
    {
        // Add the bounding box
        boundingBoxNodes_.emplace_back(node);

        // Add the bounding box's coordinates
        boundingBoxCoordinates_.insert(boundingBoxCoordinates_.end(), coordBegin, coordEnd);

        // return the index of the new node
        return boundingBoxNodes_.size() - 1;
    }

    //! Compute the bounding box of a vector of bounding boxes
    std::array<ctype, 2*dimworld>
    computeBBoxOfBBoxes_(const std::vector<ctype>& leafBoxes,
                         const std::vector<std::size_t>::iterator& begin,
                         const std::vector<std::size_t>::iterator& end)
    {
        std::array<ctype, 2*dimworld> bBoxCoords;

        // copy the iterator and get coordinates of first box
        auto it = begin;
        const auto* bFirst = leafBoxes.data() + 2*dimworld*(*it);

        for (int coordIdx = 0; coordIdx < 2*dimworld; ++coordIdx)
            bBoxCoords[coordIdx] = bFirst[coordIdx];

        // Compute min and max with the remaining boxes
        for (; it != end; ++it)
        {
            const auto* b = leafBoxes.data() + 2*dimworld*(*it);
            for (int coordIdx = 0; coordIdx < dimworld; ++coordIdx)
                if (b[coordIdx] < bBoxCoords[coordIdx])
                    bBoxCoords[coordIdx] = b[coordIdx];
            for (int coordIdx = dimworld; coordIdx < 2*dimworld; ++coordIdx)
                if (b[coordIdx] > bBoxCoords[coordIdx])
                    bBoxCoords[coordIdx] = b[coordIdx];
        }

        return bBoxCoords;
    }

    //! Compute the bounding box of a vector of bounding boxes
    std::size_t computeLongestAxis_(const std::array<ctype, 2*dimworld>& bCoords)
    {
        std::array<ctype, dimworld> axisLength;
        for (int coordIdx = 0; coordIdx < dimworld; ++coordIdx)
            axisLength[coordIdx] = bCoords[dimworld + coordIdx] - bCoords[coordIdx];

        return std::distance(axisLength.begin(), std::max_element(axisLength.begin(), axisLength.end()));
    }
};

/*!
 * \brief Check whether a point is intersectin a bounding box (dimworld == 3)
 * \param point The point
 * \param b Pointer to bounding box coordinates
 */
template<class ctype, int dimworld, typename std::enable_if_t<dimworld == 3, int> = 0>
inline bool intersectsPointBoundingBox(const Dune::FieldVector<ctype, dimworld>& point, const ctype* b)
{
    static constexpr ctype eps_ = 1.0e-7;

    using std::max;
    const auto dx = b[3] - b[0];
    const auto dy = b[4] - b[1];
    const auto dz = b[5] - b[2];
    const ctype eps = max({dx, dy, dz})*eps_;
    return (b[0] - eps <= point[0] && point[0] <= b[3] + eps &&
            b[1] - eps <= point[1] && point[1] <= b[4] + eps &&
            b[2] - eps <= point[2] && point[2] <= b[5] + eps);
}

/*!
 * \brief Check whether a point is intersectin a bounding box  (dimworld == 2)
 * \param point The point
 * \param b Pointer to bounding box coordinates
 */
template<class ctype, int dimworld, typename std::enable_if_t<dimworld == 2, int> = 0>
inline bool intersectsPointBoundingBox(const Dune::FieldVector<ctype, dimworld>& point, const ctype* b)
{
    static constexpr ctype eps_ = 1.0e-7;

    using std::max;
    const auto dx = b[2] - b[0];
    const auto dy = b[3] - b[1];
    const ctype eps = max(dx, dy)*eps_;
    return (b[0] - eps <= point[0] && point[0] <= b[2] + eps &&
            b[1] - eps <= point[1] && point[1] <= b[3] + eps);
}

/*!
 * \brief Check whether a point is intersectin a bounding box  (dimworld == 1)
 * \param point The point
 * \param b Pointer to bounding box coordinates
 */
template<class ctype, int dimworld, typename std::enable_if_t<dimworld == 1, int> = 0>
inline bool intersectsPointBoundingBox(const Dune::FieldVector<ctype, dimworld>& point, const ctype* b)
{
    static constexpr ctype eps_ = 1.0e-7;
    const ctype eps0 = eps_*(b[1] - b[0]);
    return b[0] - eps0 <= point[0] && point[0] <= b[1] + eps0;
}

/*!
 * \ingroup Geometry
 * \brief Determine if a point intersects an axis-aligned bounding box
 * The bounding box is given by the lower left corner (min) and the upper right corner (max)
 */
template<class ctype, int dimworld>
inline bool intersectsPointBoundingBox(const Dune::FieldVector<ctype, dimworld>& point,
                                       const Dune::FieldVector<ctype, dimworld>& min,
                                       const Dune::FieldVector<ctype, dimworld>& max)
{
    std::array<ctype, 2*dimworld> bBox;
    std::copy(min.begin(), min.end(), bBox.begin());
    std::copy(max.begin(), max.end(), bBox.begin()+dimworld);
    return intersectsPointBoundingBox(point, bBox.data());
}

/*!
 * \brief Check whether a bounding box is intersecting another bounding box (dimworld == 3)
 * \param a Pointer to first bounding box coordinates
 * \param b Pointer to second bounding box coordinates
 */
template<int dimworld, class ctypea, class ctypeb, typename std::enable_if_t<dimworld == 3, int> = 0>
inline bool intersectsBoundingBoxBoundingBox(const ctypea* a, const ctypeb* b)
{
    using ctype = typename Dune::PromotionTraits<ctypea, ctypeb>::PromotedType;
    static constexpr ctype eps_ = 1.0e-7;
    const ctype eps0 = eps_*std::max(b[3]-b[0], a[3]-a[0]);
    const ctype eps1 = eps_*std::max(b[4]-b[1], a[4]-a[1]);
    const ctype eps2 = eps_*std::max(b[5]-b[2], a[5]-a[2]);
    return (b[0] - eps0 <= a[3] && a[0] <= b[3] + eps0 &&
            b[1] - eps1 <= a[4] && a[1] <= b[4] + eps1 &&
            b[2] - eps2 <= a[5] && a[2] <= b[5] + eps2);

}

/*!
 * \brief Check whether a bounding box is intersecting another bounding box (dimworld == 2)
 * \param a Pointer to first bounding box coordinates
 * \param b Pointer to second bounding box coordinates
 */
template<int dimworld, class ctypea, class ctypeb, typename std::enable_if_t<dimworld == 2, int> = 0>
inline bool intersectsBoundingBoxBoundingBox(const ctypea* a, const ctypeb* b)
{
    using ctype = typename Dune::PromotionTraits<ctypea, ctypeb>::PromotedType;
    static constexpr ctype eps_ = 1.0e-7;
    const ctype eps0 = eps_*std::max(b[2]-b[0], a[2]-a[0]);
    const ctype eps1 = eps_*std::max(b[3]-b[1], a[3]-a[1]);
    return (b[0] - eps0 <= a[2] && a[0] <= b[2] + eps0 &&
            b[1] - eps1 <= a[3] && a[1] <= b[3] + eps1);
}

/*!
 * \brief Check whether a bounding box is intersecting another bounding box  (dimworld == 1)
 * \param a Pointer to first bounding box coordinates
 * \param b Pointer to second bounding box coordinates
 */
template<int dimworld, class ctypea, class ctypeb, typename std::enable_if_t<dimworld == 1, int> = 0>
inline bool intersectsBoundingBoxBoundingBox(const ctypea* a, const ctypeb* b)
{
    using ctype = typename Dune::PromotionTraits<ctypea, ctypeb>::PromotedType;
    static constexpr ctype eps_ = 1.0e-7;
    const ctype eps0 = eps_*std::max(b[1]-b[0], a[1]-a[0]);
    return b[0] - eps0 <= a[1] && a[0] <= b[1] + eps0;
}

} // end namespace Dumux

#endif
