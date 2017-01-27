// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup EmbeddedCoupling
 * \brief A class glueing two grids of different dimension geometrically
 *        Intersections are computed using axis-aligned bounding box trees
 */

#ifndef DUMUX_MULTIDIMENSION_GLUE_HH
#define DUMUX_MULTIDIMENSION_GLUE_HH

#include <iostream>
#include <fstream>
#include <string>
#include <utility>

#include <dune/common/timer.hh>
#include <dune/common/iteratorrange.hh>
#include <dune/geometry/affinegeometry.hh>

#include <dumux/common/geometrycollision.hh>
#include <dumux/common/boundingboxtree.hh>

namespace Dumux
{

namespace Properties
{
// Property forward declarations
NEW_PROP_TAG(BulkProblemTypeTag);
NEW_PROP_TAG(LowDimProblemTypeTag);
NEW_PROP_TAG(GridView);
} // namespace Properties

// forward declaration
template<class TypeTag>
class CCMultiDimensionGlue;

//! Range generator to iterate with range-based for loops over all intersections
//! as follows: for (const auto& is : intersections(glue)) { ... }
template<typename TypeTag>
Dune::IteratorRange<typename CCMultiDimensionGlue<TypeTag>::Intersections::const_iterator>
intersections(const CCMultiDimensionGlue<TypeTag>& glue)
{ return {glue.ibegin(), glue.iend()}; }

namespace Glue
{

template<class TypeTag>
class Intersection
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    // obtain the type tags of the sub problems
    using BulkProblemTypeTag = typename GET_PROP_TYPE(TypeTag, BulkProblemTypeTag);
    using LowDimProblemTypeTag = typename GET_PROP_TYPE(TypeTag, LowDimProblemTypeTag);

    using BulkGridView = typename GET_PROP_TYPE(BulkProblemTypeTag, GridView);
    using LowDimGridView = typename GET_PROP_TYPE(LowDimProblemTypeTag, GridView);

    enum {
        bulkDim = BulkGridView::dimension,
        lowDimDim = LowDimGridView::dimension,
        dimWorld = BulkGridView::dimensionworld,
        dimIs = LowDimGridView::dimension,
    };

    using BulkElement = typename BulkGridView::template Codim<0>::Entity;
    using LowDimElement = typename LowDimGridView::template Codim<0>::Entity;

    using LowDimTree = Dumux::BoundingBoxTree<LowDimGridView>;
    using BulkTree = Dumux::BoundingBoxTree<BulkGridView>;

    using Geometry = Dune::AffineGeometry<Scalar, dimIs, dimWorld>;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    Intersection(const BulkTree& bulkTree, const LowDimTree& lowDimTree)
    : bulkTree_(bulkTree), lowDimTree_(lowDimTree) {}

    //! set the intersection geometry corners
    void setCorners(const std::vector<GlobalPosition>& corners)
    {
        corners_ = corners;
        assert(corners.size() == dimIs + 1);
        type_.makeSimplex(dimIs);
    }

    //! add a pair of neighbor elements
    void addNeighbors(unsigned int bulk, unsigned int lowDim)
    {
        neighbors_[0].push_back(bulk);
        neighbors_[1].push_back(lowDim);
    }

    //! get the intersection geometry
    Geometry geometry() const
    { return Geometry(type_, corners_); }

    //! get the number of neigbors
    std::size_t neighbor(unsigned int side) const
    { return neighbors_[side].size(); }

    //! get the inside neighbor
    LowDimElement inside(unsigned int n) const
    { return lowDimTree_.entity(neighbors_[1][n]); }

    //! get the outside neighbor
    BulkElement outside(unsigned int n) const
    { return bulkTree_.entity(neighbors_[0][n]); }

private:
    std::array<std::vector<unsigned int>, 2> neighbors_;
    std::vector<GlobalPosition> corners_;
    Dune::GeometryType type_;

    const BulkTree& bulkTree_;
    const LowDimTree& lowDimTree_;
};
}

/*!
 * \ingroup EmbeddedCoupling
 * \ingroup CCModel
 * \brief Manages the coupling between bulk elements and lower dimensional elements
 *        Point sources on each integration point are computed by an AABB tree.
 *        Both domain are assumed to be discretized using a cc finite volume scheme.
 */
template<class TypeTag>
class CCMultiDimensionGlue
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);

    // obtain the type tags of the sub problems
    using BulkProblemTypeTag = typename GET_PROP_TYPE(TypeTag, BulkProblemTypeTag);
    using LowDimProblemTypeTag = typename GET_PROP_TYPE(TypeTag, LowDimProblemTypeTag);

    using BulkGridView = typename GET_PROP_TYPE(BulkProblemTypeTag, GridView);
    using LowDimGridView = typename GET_PROP_TYPE(LowDimProblemTypeTag, GridView);

    using BulkProblem = typename GET_PROP_TYPE(BulkProblemTypeTag, Problem);
    using LowDimProblem = typename GET_PROP_TYPE(LowDimProblemTypeTag, Problem);

    enum {
        bulkDim = BulkGridView::dimension,
        lowDimDim = LowDimGridView::dimension,
        dimWorld = BulkGridView::dimensionworld
    };

    enum {
        bulkIsBox = GET_PROP_VALUE(BulkProblemTypeTag, ImplicitIsBox),
        lowDimIsBox = GET_PROP_VALUE(LowDimProblemTypeTag, ImplicitIsBox)
    };

    using BulkElement = typename BulkGridView::template Codim<0>::Entity;
    using LowDimElement = typename LowDimGridView::template Codim<0>::Entity;

    using LowDimTree = Dumux::BoundingBoxTree<LowDimGridView>;
    using BulkTree = Dumux::BoundingBoxTree<BulkGridView>;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    // export intersection container type
    using Intersections = std::vector<Dumux::Glue::Intersection<TypeTag>>;

    /*!
     * \brief Constructor
     */
    CCMultiDimensionGlue(BulkProblem& bulkProblem, LowDimProblem& lowDimProblem)
    : bulkProblem_(bulkProblem), lowDimProblem_(lowDimProblem)
    {
        // Check if we are using the cellcentered method in both domains
        static_assert(!bulkIsBox && !lowDimIsBox, "Using the cell-centered glue for problems using box discretization!");
    }

    void build()
    {
        // Initialize the bulk and the lowdim bounding box tree
        const auto& bulkTree = bulkProblem().boundingBoxTree();
        const auto& lowDimTree = lowDimProblem().boundingBoxTree();

        // compute raw intersections
        Dune::Timer timer;
        auto rawIntersections = bulkTree.computeEntityCollisions(lowDimTree);
        std::cout << "Computed tree intersections in " << timer.elapsed() << std::endl;

        // create a map to check whether intersection geometries were already inserted
        std::vector<std::vector<std::vector<GlobalPosition>>> intersectionMap;
        std::vector<std::vector<unsigned int>> intersectionIndex;
        intersectionMap.resize(lowDimProblem().gridView().size(0));
        intersectionIndex.resize(lowDimProblem().gridView().size(0));

        for (const auto& rawIntersection : rawIntersections)
        {
            bool add = true;
            for (int i = 0; i < intersectionMap[rawIntersection.second()].size(); ++i)
            {
                if (rawIntersection.cornersMatch(intersectionMap[rawIntersection.second()][i]))
                {
                    add = false;
                    // only add the pair of neighbors using the insertionIndex
                    auto idx = intersectionIndex[rawIntersection.second()][i];
                    intersections_[idx].addNeighbors(rawIntersection.first(), rawIntersection.second());
                    break;
                }
            }
            if(add)
            {
                // add to the map
                intersectionMap[rawIntersection.second()].push_back(rawIntersection.corners());
                intersectionIndex[rawIntersection.second()].push_back(intersections_.size());
                // add new intersection and add the neighbors
                intersections_.emplace_back(bulkTree, lowDimTree);
                intersections_.back().setCorners(rawIntersection.corners());
                intersections_.back().addNeighbors(rawIntersection.first(), rawIntersection.second());
            }
        }
    }

    //! Return a reference to the bulk problem
    BulkProblem& bulkProblem()
    {
        return bulkProblem_;
    }

    //! Return a reference to the low dimensional problem
    LowDimProblem& lowDimProblem()
    {
        return lowDimProblem_;
    }

    //! Return begin iterator to intersection container
    typename Intersections::const_iterator ibegin() const
    { return intersections_.begin(); }

    //! Return end iterator to intersection container
    typename Intersections::const_iterator iend() const
    { return intersections_.end(); }

private:
    BulkProblem& bulkProblem_;
    LowDimProblem& lowDimProblem_;

    Intersections intersections_;
};

} // namespace Dumux

#endif
