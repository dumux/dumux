/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \brief A bounding box tree for grid point intersections
 *
 * Dumux implementation of the bounding box tree
 * adapted from implementation in FEniCS by Anders Logg
 */
#ifndef DUMUX_BOUNDINGBOXTREE_HH
#define DUMUX_BOUNDINGBOXTREE_HH

#include <memory>
#include <dune/geometry/referenceelements.hh>

namespace Dumux {

// optimized dimension-dependent methods
template <int dimworld>
class BoundingBoxTreeHelper {};

template <>
class BoundingBoxTreeHelper<3>
{
    // An epsilon for floating point operations
    static constexpr double eps_ = 1.0e-7;
public:
    // Check whether a point is in a bounding box
    static bool pointInBoundingBox(const Dune::FieldVector<double, 3>& point, unsigned int node)
    {
        const double* b = boundingBoxCoordinates_.data() + 6*node;
        const double eps0 = eps_*(b[3] - b[0]);
        const double eps1 = eps_*(b[4] - b[1]);
        const double eps2 = eps_*(b[5] - b[2]);
        return (b[0] - eps0 <= point[0] && point[0] <= b[3] + eps0 &&
        	    b[1] - eps1 <= point[1] && point[1] <= b[4] + eps1 &&
        	    b[2] - eps2 <= point[2] && point[2] <= b[5] + eps2);
    }

    // Compute the bounding box of a vector of bounding boxes
    static void computeBBoxOfBBoxes(double* bBox,
                                    std::size_t& axis,
                                    const std::vector<double>& leafBoxes,
                                    const std::vector<unsigned int>::iterator& begin,
                                    const std::vector<unsigned int>::iterator& end)
    {
        // Copy the iterator and get coordinates of first box
        auto it = begin;
        const double* b = leafBoxes.data() + 6*(*it);
        // Maybe write out the loop for optimization
        for (std::size_t coordIdx = 0; coordIdx < 6; ++coordIdx)
            bBox[coordIdx] = b[coordIdx];

        // Compute min and max with the remaining boxes
        for (; it != end; ++it)
        {
            const double* b = leafBoxes.data() + 6*(*it);
            // Maybe write out the loop for optimization
            for (std::size_t coordIdx = 0; coordIdx < 6; ++coordIdx)
                if (b[coordIdx] < bBox[coordIdx]) bBox[coordIdx] = b[coordIdx];
        }

        // Compute the longest axis
        const double x = bBox[3] - bBox[0];
        const double y = bBox[4] - bBox[1];
        const double z = bBox[5] - bBox[2];

        if (x > y && x > z)
            axis = 0;
        else if (y > z)
            axis = 1;
        else
            axis = 2;
    }

    // Sort the bounding boxes along the longest axis
    static void sortBoundingBoxes(std::size_t axis,
                                  const std::vector<double>& leafBoxes,
                                  const std::vector<unsigned int>::iterator& begin,
                                  const std::vector<unsigned int>::iterator& middle,
                                  const std::vector<unsigned int>::iterator& end)
    {
        switch (axis)
        {
            case 0:
                std::nth_element(begin, middle, end, lessXBox(leafBoxes));
            case 1:
                std::nth_element(begin, middle, end, lessYBox(leafBoxes));
            default:
                std::nth_element(begin, middle, end, lessZBox(leafBoxes));
        }
    }

    // Comparison operators for sorting bounding boxes. There are sorted by their
    // mid points along the longest axis
    struct lessXBox
    {
        const std::vector<double>& bBoxes;
        lessXBox(const std::vector<double>& bBoxes): bBoxes(bBoxes) {}
        inline bool operator()(unsigned int i, unsigned int j)
        {
            const double* bi = bBoxes.data() + 6*i;
            const double* bj = bBoxes.data() + 6*j;
            return bi[0] + bi[3] < bj[0] + bj[3];
        }
    };

    struct lessYBox
    {
        const std::vector<double>& bBoxes;
        lessYBox(const std::vector<double>& bBoxes): bBoxes(bBoxes) {}
        inline bool operator()(unsigned int i, unsigned int j)
        {
            const double* bi = bBoxes.data() + 6*i;
            const double* bj = bBoxes.data() + 6*j;
            return bi[1] + bi[4] < bj[1] + bj[4];
        }
    };

    struct lessZBox
    {
        const std::vector<double>& bBoxes;
        lessZBox(const std::vector<double>& bBoxes): bBoxes(bBoxes) {}
        inline bool operator()(unsigned int i, unsigned int j)
        {
            const double* bi = bBoxes.data() + 6*i;
            const double* bj = bBoxes.data() + 6*j;
            return bi[2] + bi[5] < bj[2] + bj[5];
        }
    };
};

template <>
class BoundingBoxTreeHelper<2>
{
    // An epsilon for floating point operations
    static constexpr double eps_ = 1.0e-7;
public:
    // Check whether a point is in a bounding box
    static bool pointInBoundingBox(const Dune::FieldVector<double, 2>& point, unsigned int node)
    {
        const double* b = boundingBoxCoordinates_.data() + 4*node;
        const double eps0 = eps_*(b[3] - b[0]);
        const double eps1 = eps_*(b[4] - b[1]);
        return (b[0] - eps0 <= point[0] && point[0] <= b[2] + eps0 &&
        	    b[1] - eps1 <= point[1] && point[1] <= b[3] + eps1);
    }

    // Compute the bounding box of a vector of bounding boxes
    static void computeBBoxOfBBoxes(double* bBox,
                                    std::size_t& axis,
                                    const std::vector<double>& leafBoxes,
                                    const std::vector<unsigned int>::iterator& begin,
                                    const std::vector<unsigned int>::iterator& end)
    {
        // Copy the iterator and get coordinates of first box
        auto it = begin;
        const double* b = leafBoxes.data() + 4*(*it);
        // Maybe write out the loop for optimization
        for (std::size_t coordIdx = 0; coordIdx < 4; ++coordIdx)
            bBox[coordIdx] = b[coordIdx];

        // Compute min and max with the remaining boxes
        for (; it != end; ++it)
        {
            const double* b = leafBoxes.data() + 4*(*it);
            // Maybe write out the loop for optimization
            for (std::size_t coordIdx = 0; coordIdx < 4; ++coordIdx)
                if (b[coordIdx] < bBox[coordIdx]) bBox[coordIdx] = b[coordIdx];
        }

        // Compute the longest axis
        const double x = bBox[2] - bBox[0];
        const double y = bBox[3] - bBox[1];

        if (x > y)
            axis = 0;
        else
            axis = 1;
    }

    // Sort the bounding boxes along the longest axis
    static void sortBoundingBoxes(std::size_t axis,
                                  const std::vector<double>& leafBoxes,
                                  const std::vector<unsigned int>::iterator& begin,
                                  const std::vector<unsigned int>::iterator& middle,
                                  const std::vector<unsigned int>::iterator& end)
    {
        if (axis == 0)
            std::nth_element(begin, middle, end, lessXBox(leafBoxes));
        else
            std::nth_element(begin, middle, end, lessYBox(leafBoxes));
    }

    // Comparison operators for sorting bounding boxes. There are sorted by their
    // mid points along the longest axis
    struct lessXBox
    {
        const std::vector<double>& bBoxes;
        lessXBox(const std::vector<double>& bBoxes): bBoxes(bBoxes) {}
        inline bool operator()(unsigned int i, unsigned int j)
        {
            const double* bi = bBoxes.data() + 4*i;
            const double* bj = bBoxes.data() + 4*j;
            return bi[0] + bi[2] < bj[0] + bj[2];
        }
    };

    struct lessYBox
    {
        const std::vector<double>& bBoxes;
        lessYBox(const std::vector<double>& bBoxes): bBoxes(bBoxes) {}
        inline bool operator()(unsigned int i, unsigned int j)
        {
            const double* bi = bBoxes.data() + 4*i;
            const double* bj = bBoxes.data() + 4*j;
            return bi[1] + bi[3] < bj[1] + bj[3];
        }
    };
};

//! An index to element map
template <class GridView>
class IndexToElementMap
  : public std::vector<typename GridView::Traits::Grid::template Codim<0>::EntitySeed>
{
    typedef typename GridView::template Codim<0>::Entity Element;
public:
    IndexToElementMap(const GridView& gridView)
      : grid(gridView.grid()) {}

    template<class T>
    Element entity(T&& t)
    { return grid_.entity(*this[std::forward<T>(t)]); }

private:
    const GridView::Traits::Grid& grid_;
};

//! The bounding box class. Implements an axis-aligned bounding box tree for grids.
template <class GridView>
class BoundingBoxTree
{
    static const int dim = GridView::dimension;
    static const int dimworld = GridView::dimensionworld;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename Dune::ReferenceElements<double, dim> ReferenceElements;
    typedef typename Dune::ReferenceElement<double, dim> ReferenceElement;

public:
    // Default Constructor
    BoundingBoxTree() {}

    // Constructor with gridView
    BoundingBoxTree(const GridView& leafGridView)
    { build(leafGridView); }

    // Destructor
    ~BoundingBoxTree() {}

    // Build up bounding box tree for a grid
    void build(const GridView& leafGridView)
    {
        // clear data if any
        clear_();

        // Create bounding boxes for all elements
        const unsigned int numLeaves = leafGridView.size(0);
        std::vector<double> leafBoxes(2*dimworld*numLeaves); // min and max for all dims
        indexToElementMap_ = std::make_shared<IndexToElementMap<GridView> >(leafGridView);
        indexToElementMap_->resize(numLeaves);
        for (auto&& element : elements(leafGridView))
        {
            unsigned int eIdx = leafGridView.indexSet().index(element);
            computeEntityBoundingBox_(leafBoxes.data() + 2*dimworld*eIdx, element);
            indexToElementMap_->push_back(element.seed());
        }

        // Create the leaf partition (to be sorted)
        std::vector<unsigned int> leafPartition(numLeaves);
        for (unsigned int i = 0; i < numLeaves; ++i)
            leafPartition[i] = i;

        // Recursively build the bounding box tree starting from the leaves
        build_(leafBoxes, leafPartition.begin(), leafPartition.end());

        // We are done, log output
        std::cout << "Computed bounding box tree with " << numBoundingBoxes_()
                  << " nodes for " << numLeaves << " grid entites." << std::endl;
    }

    // Compute all intersections between entities and a point
    std::vector<unsigned int> computeEntityCollisions(const Dune::FieldVector<double, dimworld>& point) const
    {
        // Call the recursive find function to find candidates
        std::vector<unsigned int> entities;
        computeCollisions_(point, numBoundingBoxes_() - 1, entities);
        return entities;
    }

private:

    // Bounding box data. Leaf nodes are indicated by setting child_0 to
    // the node itself and child_1 is the index of the entity in the bounding box.
    struct BoundingBox
    {
        unsigned int child_0;
        unsigned int child_1;
    };

    // Vector of bounding boxes
    std::vector<BoundingBox> boundingBoxes_;

    // Vector of bounding box coordinates
    std::vector<double> boundingBoxCoordinates_;

    // Shared pointer to the index to element map
    std::shared_ptr<IndexToElementMap<GridView> > indexToElementMap_;

    // Clear all data
    void clear_()
    {
        boundingBoxes_.clear();
        boundingBoxCoordinates_.clear();
        if(indexToElementMap_) indexToElementMap_->clear();
    }

    // Build bounding box tree for all entities recursively
    unsigned int build_(const std::vector<double>& leafBoxes,
                        const std::vector<unsigned int>::iterator& begin,
                        const std::vector<unsigned int>::iterator& end)
    {
        assert(begin > end);

        // Create empty bounding box
        BoundingBox bBox;

        // If we reached the leaf
        if (end - begin == 1)
        {
            // Get the bounding box coordinates for the leaf
            const unsigned int eIdx = *begin;
            const double* b = leafBoxes.data() + 2*dimworld*eIdx;

            // Store the data in the bounding box
            bBox.child_0 = numBoundingBoxes_();
            bBox.child_1 = eIdx;
            return addBoundingBox_(bBox, b);
        }

        // Compute the bounding box of all bounding boxes in the range [begin, end]
        double b[dimworld];
        std::size_t axis;
        BoundingBoxTreeHelper<dimworld>::computeBBoxOfBBoxes(b, axis, leafBoxes, begin, end);

        // Sort bounding boxes along the longest axis
        std::vector<unsigned int>::iterator middle = begin + (end - begin) / 2;
        BoundingBoxTreeHelper<dimworld>::sortBoundingBoxes(axis, leafBoxes, begin, middle, end);

        // Split the bounding boxes into two branches and call build recursively
        bBox.child_0 = build_(leafBoxes, begin, middle);
        bBox.child_1 = build_(leafBoxes, middle, end);

        // Store the bounding box data. The root will be added at the end.
        return addBoundingBox_(bBox, b);
    }

    // Compute collisions with point recursively
    void computeCollisions_(const Dune::FieldVector<double, dimworld>& point,
                            unsigned int node,
                            std::vector<unsigned int>& entities)
    {
        // Get the bounding box for the current node
        const BoundingBox& bBox = getBoundingBox_(node);

        // if the point is not in the bounding box we can stop
        if (!BoundingBoxTreeHelper<dimworld>::pointInBoundingBox(point, node))
            return;

        // We know now it's inside. If the box is a leaf add it.
        else if (isLeaf_(bBox, node))
        {
            const unsigned int eIdx = bBox.child_1;
            auto geometry = (indexToElementMap_->entity(eIdx)).geometry();
            const ReferenceElement &refElement = ReferenceElements::general(geometry.type());
            if (refElement.checkInside(geometry.local(point)))
                entities.push_back(eIdx);
        }

        // No leaf. Check both children.
        else
        {
            computeCollisions_(point, bBox.child_0, entities);
            computeCollisions_(point, bBox.child_1, entities);
        }
    }

    // Add a new bounding box to the tree
    inline unsigned int addBoundingBox_(const BoundingBox& bBox,
                                        const double* b)
    {
        // Add the bounding box
        boundingBoxes_.push_back(bBox);

        // Add the bounding box's coordinates
        for (std::size_t i = 0; i < 2*dimworld; ++i)
            boundingBoxCoordinates_.push_back(b[i]);

        // return the index of the new node
        return boundingBoxes_.size() - 1;
    }

    // Get an existing bounding box for a given node
    inline const BoundingBox& getBoundingBox_(unsigned int node) const
    { return boundingBoxes_[node]; }

    // Get the number of bounding boxes currently in the tree
    inline std::size_t numBoundingBoxes_() const
    { return boundingBoxes_.size(); }

    // Check whether a bounding box is a leaf node
    // Leaf nodes have itself as child_0
    inline bool isLeaf_(const BoundingBox& bBox, unsigned int node) const
    { return bBox.child_0 == node; }

    // Compute the bounding box of a grid entity
    void computeEntityBoundingBox_(double* b, const Entity& entity) const
    {
        // get the bounding box coordinates
        double* xMin = b;
        double* xMax = b + dimworld;

        // get mesh entity data
        auto geometry = entity.geometry();

        // Get coordinates of first vertex
        auto corner = geometry.corner(0);
        for (std::size_t dimIdx = 0; dimIdx < dimworld; ++dimIdx)
            xMin[dimIdx] = xMax[dimIdx] = corner[dimIdx];

        // Compute the min and max over the remaining vertices
        for (std::size_t vLocalIdx = 1; vLocalIdx < entity.subEntities(dim); ++vLocalIdx)
        {
            corner = geometry.corner(vLocalIdx);
            for (std::size_t dimIdx = 0; dimIdx < dimworld; ++dimIdx)
            {
                xMin[dimIdx] = std::min(xMin[dimIdx], corner(vLocalIdx)[dimIdx]);
                xMax[dimIdx] = std::max(xMax[dimIdx], corner(vLocalIdx)[dimIdx]);
            }
        }
    }
};

} // end namespace Dumux

#endif
