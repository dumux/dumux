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
 * \ingroup Geometry
 * \brief A class representing the intersection entites two geometric entity sets
 */

#ifndef DUMUX_GEOMETRY_INTERSECTION_ENTITY_SET_HH
#define DUMUX_GEOMETRY_INTERSECTION_ENTITY_SET_HH

#include <algorithm>
#include <cassert>
#include <iostream>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/indices.hh>
#include <dune/common/timer.hh>
#include <dune/common/iteratorrange.hh>
#include <dune/common/promotiontraits.hh>
#include <dune/common/reservedvector.hh>
#include <dune/geometry/affinegeometry.hh>
#include <dune/geometry/type.hh>

#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/geometry/intersectingentities.hh>

namespace Dumux {

/*!
 * \ingroup Geometry
 * \brief A class representing the intersection entites two geometric entity sets
 */
template<class DomainEntitySet, class TargetEntitySet>
class IntersectionEntitySet
{
    using DomainTree = BoundingBoxTree<DomainEntitySet>;
    using TargetTree = BoundingBoxTree<TargetEntitySet>;

    using ctype = typename Dune::PromotionTraits<typename DomainEntitySet::ctype, typename TargetEntitySet::ctype>::PromotedType;

    static constexpr int dimWorld = DomainEntitySet::dimensionworld;
    static_assert(dimWorld == int(TargetEntitySet::dimensionworld), "Entity sets must have the same world dimension");

    using GlobalPosition = Dune::FieldVector<ctype, dimWorld>;

    static constexpr int dimDomain = DomainEntitySet::Entity::Geometry::mydimension;
    static constexpr int dimTarget = TargetEntitySet::Entity::Geometry::mydimension;
    static constexpr bool isMixedDimensional = dimDomain != dimTarget;

    /*!
     * \brief A class representing an intersection entity
     */
    class IntersectionEntity
    {
        static constexpr int dimIs = std::min(dimDomain, dimTarget);
        using Geometry = Dune::AffineGeometry<ctype, dimIs, dimWorld>; // geometries are always simplices

        // we can only have multiple neighbors in the mixeddimensional case and then only for the side with the largest dimension
        using IndexStorage = std::pair<std::conditional_t<dimDomain <= dimTarget, Dune::ReservedVector<std::size_t, 1>, std::vector<std::size_t>>,
                                       std::conditional_t<dimTarget <= dimDomain, Dune::ReservedVector<std::size_t, 1>, std::vector<std::size_t>>>;

        static constexpr auto domainIdx = Dune::index_constant<0>{};
        static constexpr auto targetIdx = Dune::index_constant<1>{};

    public:
        IntersectionEntity(const DomainTree& domainTree, const TargetTree& targetTree)
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
            if (numDomainNeighbors() == 0 && numTargetNeighbors() == 0)
            {
                std::get<domainIdx>(neighbors_).push_back(domain);
                std::get<targetIdx>(neighbors_).push_back(target);
            }
            else if (dimDomain > dimTarget)
                std::get<domainIdx>(neighbors_).push_back(domain);

            else if (dimTarget > dimDomain)
                std::get<targetIdx>(neighbors_).push_back(target);

            else
                DUNE_THROW(Dune::InvalidStateException, "Cannot add more than one neighbor per side for equidimensional intersection!");
        }

        //! get the intersection geometry
        Geometry geometry() const
        { return Geometry(Dune::GeometryTypes::simplex(dimIs), corners_); }

        //! get the number of domain neighbors of this intersection
        std::size_t numDomainNeighbors() const
        { return std::get<domainIdx>(neighbors_).size(); }

        //! get the number of target neighbors of this intersection
        std::size_t numTargetNeighbors() const
        { return std::get<targetIdx>(neighbors_).size(); }

        //! get the nth domain neighbor entity
        typename DomainEntitySet::Entity domainEntity(unsigned int n = 0) const
        { return domainTree_.entitySet().entity(std::get<domainIdx>(neighbors_)[n]); }

        //! get the nth target neighbor entity
        typename TargetEntitySet::Entity targetEntity(unsigned int n = 0) const
        { return targetTree_.entitySet().entity(std::get<targetIdx>(neighbors_)[n]); }

    private:
        IndexStorage neighbors_;
        std::vector<GlobalPosition> corners_;

        const DomainTree& domainTree_;
        const TargetTree& targetTree_;
    };

    using Intersections = std::vector<IntersectionEntity>;

public:
    //! make intersection entity type available
    using Entity = IntersectionEntity;
    //! make entity iterator type available
    using EntityIterator = typename Intersections::const_iterator;

    /*!
     * \brief Default constructor
     */
    IntersectionEntitySet() = default;

    /*!
     * \brief Build intersections
     */
    void build(std::shared_ptr<const DomainEntitySet> domainSet, std::shared_ptr<const TargetEntitySet> targetSet)
    {
        domainTree_ = std::make_shared<DomainTree>(domainSet);
        targetTree_ = std::make_shared<TargetTree>(targetSet);
        build(*domainTree_, *targetTree_);
    }

    /*!
     * \brief Build intersections
     */
    void build(std::shared_ptr<const DomainTree> domainTree, std::shared_ptr<const TargetTree> targetTree)
    {
        // make sure the tree don't get out of scope
        domainTree_ = domainTree;
        targetTree_ = targetTree;
        build(*domainTree_, *targetTree_);
    }

    /*!
     * \brief Build intersections
     * \note If you call this, make sure the bounding box tree stays alive for the life-time of this object
     */
    void build(const DomainTree& domainTree, const TargetTree& targetTree)
    {
        Dune::Timer timer;

        // compute raw intersections
        const auto rawIntersections = intersectingEntities(domainTree, targetTree);

        // create a map to check whether intersection geometries were already inserted
        // Note that this can only occur if the grids have different dimensionality.
        // If this is the case, we keep track of the intersections using the indices of the lower-
        // dimensional entities which is identical for all geometrically identical intersections.
        std::vector<std::vector<std::vector<GlobalPosition>>> intersectionMap;
        std::vector<std::vector<std::size_t>> intersectionIndex;
        if constexpr (isMixedDimensional)
        {
            const auto numLowDimEntities = dimTarget < dimDomain ? targetTree.entitySet().size()
                                                                 : domainTree.entitySet().size();
            intersectionMap.resize(numLowDimEntities);
            intersectionIndex.resize(numLowDimEntities);
        }

        // reserve memory for storing the intersections. In case of grids of
        // different dimensionality this might be an overestimate. We get rid
        // of the overhead memory at the end of this function though.
        intersections_.clear();
        intersections_.reserve(rawIntersections.size());

        for (const auto& rawIntersection : rawIntersections)
        {
            bool add = true;

            // Check if intersection was already inserted.
            // In this case we only add new neighbor information as the geometry is identical.
            if constexpr (isMixedDimensional)
            {
                const auto lowDimNeighborIdx = getLowDimNeighborIdx_(rawIntersection);
                for (int i = 0; i < intersectionMap[lowDimNeighborIdx].size(); ++i)
                {
                    if (rawIntersection.cornersMatch(intersectionMap[lowDimNeighborIdx][i]))
                    {
                        add = false;
                        // only add the pair of neighbors using the insertionIndex
                        auto idx = intersectionIndex[lowDimNeighborIdx][i];
                        intersections_[idx].addNeighbors(rawIntersection.first(), rawIntersection.second());
                        break;
                    }
                }
            }

            if(add)
            {
                // maybe add to the map
                if constexpr (isMixedDimensional)
                {
                    intersectionMap[getLowDimNeighborIdx_(rawIntersection)].push_back(rawIntersection.corners());
                    intersectionIndex[getLowDimNeighborIdx_(rawIntersection)].push_back(intersections_.size());
                }

                // add new intersection and add the neighbors
                intersections_.emplace_back(domainTree, targetTree);
                intersections_.back().setCorners(rawIntersection.corners());
                intersections_.back().addNeighbors(rawIntersection.first(), rawIntersection.second());
            }
        }

        intersections_.shrink_to_fit();
        std::cout << "Computed " << size() << " intersection entities in " << timer.elapsed() << std::endl;
    }

    //! return begin iterator to intersection container
    typename Intersections::const_iterator ibegin() const
    { return intersections_.begin(); }

    //! return end iterator to intersection container
    typename Intersections::const_iterator iend() const
    { return intersections_.end(); }

    //! the number of intersections
    std::size_t size() const
    { return intersections_.size(); }

    /*!
     * \brief Range generator to iterate with range-based for loops over all intersections
     *        as follows: for (const auto& is : intersections(glue)) { ... }
     * \note free function
     */
    friend Dune::IteratorRange<EntityIterator> intersections(const IntersectionEntitySet& set)
    { return {set.ibegin(), set.iend()}; }

private:
    template<class RawIntersection,
             bool enable = isMixedDimensional, std::enable_if_t<enable, int> = 0>
    auto getLowDimNeighborIdx_(const RawIntersection& is)
    {
        if constexpr (dimTarget < dimDomain)
            return is.second();
        else
            return is.first();
    }

    Intersections intersections_;

    std::shared_ptr<const DomainTree> domainTree_;
    std::shared_ptr<const TargetTree> targetTree_;
};

} // end namespace Dumux

#endif
