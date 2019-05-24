// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup MultiDomain
 * \brief A class glueing two grids of potentially different dimension geometrically.
 *        Intersections are computed using axis-aligned bounding box trees
 */

#ifndef DUMUX_MULTIDOMAIN_GLUE_HH
#define DUMUX_MULTIDOMAIN_GLUE_HH

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
template<class DomainGridView, class TargetGridView, class DomainMapper, class TargetMapper>
class MultiDomainGlue;

/*!
 * \ingroup MultiDomain
 * \brief Range generator to iterate with range-based for loops over all intersections
 *        as follows: for (const auto& is : intersections(glue)) { ... }
 */
template<class DomainGridView, class TargetGridView, class DomainMapper, class TargetMapper>
Dune::IteratorRange<typename MultiDomainGlue<DomainGridView, TargetGridView, DomainMapper, TargetMapper>::Intersections::const_iterator>
intersections(const MultiDomainGlue<DomainGridView, TargetGridView, DomainMapper, TargetMapper>& glue)
{ return {glue.ibegin(), glue.iend()}; }

namespace Glue {

/*!
 * \ingroup MultiDomain
 * \brief An intersection object representing the intersection
 *        between two grid entites of different grids
 */
template<class DomainGridView, class TargetGridView, class DomainMapper, class TargetMapper>
class Intersection
{
    using DomainElement = typename DomainGridView::template Codim<0>::Entity;
    using TargetElement = typename TargetGridView::template Codim<0>::Entity;

    using DomainEntitySet = GridViewGeometricEntitySet<DomainGridView, 0, DomainMapper>;
    using DomainTree = BoundingBoxTree<DomainEntitySet>;

    using TargetEntitySet = GridViewGeometricEntitySet<TargetGridView, 0, TargetMapper>;
    using TargetTree = BoundingBoxTree<TargetEntitySet>;

    static constexpr int dimWorld = DomainGridView::dimensionworld;
    static_assert(dimWorld == int(TargetGridView::dimensionworld), "Grids must have the same world dimension");

    static constexpr int dimDomain = DomainGridView::dimension;
    static constexpr int dimTarget = TargetGridView::dimension;
    static constexpr int dimIs = std::min(dimDomain, dimTarget);

    using Scalar = typename DomainGridView::ctype;
    using Geometry = Dune::AffineGeometry<Scalar, dimIs, dimWorld>;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    Intersection(const DomainTree& domainTree, const TargetTree& targetTree)
    : domainTree_(domainTree)
    , targetTree_(targetTree)
    {}

    //! set the intersection geometry corners
    void setCorners(const std::vector<GlobalPosition>& corners)
    {
        corners_ = corners;
        assert(corners.size() == dimIs + 1); // Only simplex intersections are allowed
    }

    //! add a pair of neighbor elements
    void addNeighbors(std::size_t domain, std::size_t target)
    {
        neighbors_[0].push_back(domain);
        neighbors_[1].push_back(target);
    }

    //! get the intersection geometry
    Geometry geometry() const
    { return Geometry(Dune::GeometryTypes::simplex(dimIs), corners_); }

    //! get the number of neigbors
    std::size_t neighbor(unsigned int side) const
    { return neighbors_[side].size(); }

    //! get the inside (domain) neighbor
    DomainElement inside(unsigned int n) const
    { return domainTree_.entitySet().entity(neighbors_[0][n]); }

    //! get the outside (target) neighbor
    TargetElement outside(unsigned int n) const
    { return targetTree_.entitySet().entity(neighbors_[1][n]); }

private:
    std::array<std::vector<std::size_t>, 2> neighbors_;
    std::vector<GlobalPosition> corners_;

    const DomainTree& domainTree_;
    const TargetTree& targetTree_;
};

} // end namespace Glue

/*!
 * \ingroup MultiDomain
 * \brief Manages the coupling between domain elements and lower dimensional elements
 *        Point sources on each integration point are computed by an AABB tree.
 *        Both domain are assumed to be discretized using a cc finite volume scheme.
 */
template<class DomainGridView, class TargetGridView,
         class DomainMapper = Dune::MultipleCodimMultipleGeomTypeMapper<DomainGridView>,
         class TargetMapper = Dune::MultipleCodimMultipleGeomTypeMapper<TargetGridView>>
class MultiDomainGlue
{
    using DomainEntitySet = GridViewGeometricEntitySet<DomainGridView, 0, DomainMapper>;
    using DomainTree = BoundingBoxTree<DomainEntitySet>;

    using TargetEntitySet = GridViewGeometricEntitySet<TargetGridView, 0, TargetMapper>;
    using TargetTree = BoundingBoxTree<TargetEntitySet>;

    using Scalar = typename DomainGridView::ctype;
    enum { dimWorld = DomainGridView::dimensionworld };
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    // export intersection container type
    using Intersections = std::vector<Glue::Intersection<DomainGridView, TargetGridView, DomainMapper, TargetMapper>>;


    /*!
     * \brief Default constructor
     */
    MultiDomainGlue() = default;

    /*!
     * \brief Constructor
     */
    MultiDomainGlue(const DomainTree& domainTree, const TargetTree& targetTree)
    { build(domainTree, targetTree); }

    void build(const DomainTree& domainTree, const TargetTree& targetTree)
    {
        Dune::Timer timer;

        // compute raw intersections
        auto rawIntersections = intersectingEntities(domainTree, targetTree);

        // create a map to check whether intersection geometries were already inserted
        std::vector<std::vector<std::vector<GlobalPosition>>> intersectionMap;
        std::vector<std::vector<std::size_t>> intersectionIndex;
        intersectionMap.resize(targetTree.entitySet().size());
        intersectionIndex.resize(targetTree.entitySet().size());
        intersections_.clear();

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
                intersections_.emplace_back(domainTree, targetTree);
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

/*!
 * \ingroup MultiDomain
 * \brief Creates the glue object containing the intersections
 *        between two grids obtained from given grid geometries.
 * \param domainGridGeometry The grid geometry of the domain
 * \param targetGridGeometry The grid geometry of the target domain
 * \return The glue object containing the intersections
 */
template<class DomainGG, class TargetGG>
MultiDomainGlue< typename DomainGG::GridView, typename TargetGG::GridView,
                 typename DomainGG::ElementMapper, typename TargetGG::ElementMapper >
makeGlue(const DomainGG& domainGridGeometry, const TargetGG& targetGridGeometry)
{
    return {domainGridGeometry.boundingBoxTree(), targetGridGeometry.boundingBoxTree()};
}

} // end namespace Dumux

#endif
