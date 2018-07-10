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
 * \ingroup MultiDomain
 * \ingroup MixedDimension
 * \brief A class glueing two grids of different dimension geometrically
 *        Intersections are computed using axis-aligned bounding box trees
 */

#ifndef DUMUX_MULTIDOMAIN_MIXEDDIMENSION_GLUE_HH
#define DUMUX_MULTIDOMAIN_MIXEDDIMENSION_GLUE_HH

#include <iostream>
#include <fstream>
#include <string>
#include <utility>

#include <dune/common/timer.hh>
#include <dune/common/iteratorrange.hh>
#include <dune/geometry/affinegeometry.hh>
#include <dune/geometry/type.hh>

#include <dumux/common/geometry/geometricentityset.hh>
#include <dumux/common/geometry/boundingboxtree.hh>
#include <dumux/common/geometry/intersectingentities.hh>

namespace Dumux {

// forward declaration
template<class BulkGridView, class LowDimGridView, class BulkMapper, class LowDimMapper>
class MixedDimensionGlue;

/*!
 * \ingroup MultiDomain
 * \ingroup MixedDimension
 * \brief Range generator to iterate with range-based for loops over all intersections
 *        as follows: for (const auto& is : intersections(glue)) { ... }
 */
template<class BulkGridView, class LowDimGridView, class BulkMapper, class LowDimMapper>
Dune::IteratorRange<typename MixedDimensionGlue<BulkGridView, LowDimGridView, BulkMapper, LowDimMapper>::Intersections::const_iterator>
intersections(const MixedDimensionGlue<BulkGridView, LowDimGridView, BulkMapper, LowDimMapper>& glue)
{ return {glue.ibegin(), glue.iend()}; }

namespace Glue {

/*!
 * \ingroup MultiDomain
 * \ingroup MixedDimension
 * \brief An intersection object representing an intersection
 *        between two grid entites of different grids
 */
template<class BulkGridView, class LowDimGridView, class BulkMapper, class LowDimMapper>
class Intersection
{
    using BulkElement = typename BulkGridView::template Codim<0>::Entity;
    using LowDimElement = typename LowDimGridView::template Codim<0>::Entity;

    using BulkEntitySet = GridViewGeometricEntitySet<BulkGridView, 0, BulkMapper>;
    using BulkTree = BoundingBoxTree<BulkEntitySet>;

    using LowDimEntitySet = GridViewGeometricEntitySet<LowDimGridView, 0, LowDimMapper>;
    using LowDimTree = BoundingBoxTree<LowDimEntitySet>;

    enum {
        dimWorld = BulkGridView::dimensionworld,
        dimIs = LowDimGridView::dimension,
    };

    using Scalar = typename BulkGridView::ctype;
    using Geometry = Dune::AffineGeometry<Scalar, dimIs, dimWorld>;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    Intersection(const BulkTree& bulkTree, const LowDimTree& lowDimTree)
    : bulkTree_(bulkTree)
    , lowDimTree_(lowDimTree)
    {}

    //! set the intersection geometry corners
    void setCorners(const std::vector<GlobalPosition>& corners)
    {
        corners_ = corners;
        assert(corners.size() == dimIs + 1); // Only simplex intersections are allowed
    }

    //! add a pair of neighbor elements
    void addNeighbors(std::size_t bulk, std::size_t lowDim)
    {
        neighbors_[0].push_back(bulk);
        neighbors_[1].push_back(lowDim);
    }

    //! get the intersection geometry
    Geometry geometry() const
    { return Geometry(Dune::GeometryTypes::simplex(dimIs), corners_); }

    //! get the number of neigbors
    std::size_t neighbor(unsigned int side) const
    { return neighbors_[side].size(); }

    //! get the inside neighbor
    LowDimElement inside(unsigned int n) const
    { return lowDimTree_.entitySet().entity(neighbors_[1][n]); }

    //! get the outside neighbor
    BulkElement outside(unsigned int n) const
    { return bulkTree_.entitySet().entity(neighbors_[0][n]); }

private:
    std::array<std::vector<std::size_t>, 2> neighbors_;
    std::vector<GlobalPosition> corners_;

    const BulkTree& bulkTree_;
    const LowDimTree& lowDimTree_;
};

} // end namespace Glue

/*!
 * \ingroup MixedDimension
 * \brief Manages the coupling between bulk elements and lower dimensional elements
 *        Point sources on each integration point are computed by an AABB tree.
 *        Both domain are assumed to be discretized using a cc finite volume scheme.
 */
template<class BulkGridView, class LowDimGridView,
         class BulkMapper = Dune::MultipleCodimMultipleGeomTypeMapper<BulkGridView>,
         class LowDimMapper = Dune::MultipleCodimMultipleGeomTypeMapper<LowDimGridView>>
class MixedDimensionGlue
{
    using BulkEntitySet = GridViewGeometricEntitySet<BulkGridView, 0, BulkMapper>;
    using BulkTree = BoundingBoxTree<BulkEntitySet>;

    using LowDimEntitySet = GridViewGeometricEntitySet<LowDimGridView, 0, LowDimMapper>;
    using LowDimTree = BoundingBoxTree<LowDimEntitySet>;

    using Scalar = typename BulkGridView::ctype;
    enum { dimWorld = BulkGridView::dimensionworld };
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    // export intersection container type
    using Intersections = std::vector<Glue::Intersection<BulkGridView, LowDimGridView, BulkMapper, LowDimMapper>>;


    /*!
     * \brief Default constructor
     */
    MixedDimensionGlue() = default;

    /*!
     * \brief Constructor
     */
    MixedDimensionGlue(const BulkTree& bulkTree, const LowDimTree& lowDimTree)
    { build(bulkTree, lowDimTree); }

    void build(const BulkTree& bulkTree, const LowDimTree& lowDimTree)
    {
        Dune::Timer timer;

        // compute raw intersections
        auto rawIntersections = intersectingEntities(bulkTree, lowDimTree);

        // create a map to check whether intersection geometries were already inserted
        std::vector<std::vector<std::vector<GlobalPosition>>> intersectionMap;
        std::vector<std::vector<std::size_t>> intersectionIndex;
        intersectionMap.resize(lowDimTree.entitySet().size());
        intersectionIndex.resize(lowDimTree.entitySet().size());

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

        std::cout << "Computed tree intersections in " << timer.elapsed() << std::endl;
    }

    //! Return begin iterator to intersection container
    typename Intersections::const_iterator ibegin() const
    { return intersections_.begin(); }

    //! Return end iterator to intersection container
    typename Intersections::const_iterator iend() const
    { return intersections_.end(); }

    //! the number of intersections
    std::size_t size() const
    { return intersections_.size(); }

private:
    Intersections intersections_;
};

} // end namespace Dumux

#endif
