// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Geometry
 * \brief A class representing the intersection entities of two distributed geometric entity sets
 */
#ifndef DUMUX_GEOMETRY_DISTRIBUTED_INTERSECTION_ENTITY_SET_HH
#define DUMUX_GEOMETRY_DISTRIBUTED_INTERSECTION_ENTITY_SET_HH

#include <algorithm>
#include <memory>
#include <utility>
#include <vector>

#include <dune/common/timer.hh>
#include <dune/common/iteratorrange.hh>
#include <dune/common/iteratorfacades.hh>
#include <dune/geometry/affinegeometry.hh>
#include <dune/geometry/type.hh>

#include <dumux/geometry/distributedboundingboxtree.hh>
#include <dumux/geometry/distributedintersectingentities.hh>

namespace Dumux {

/*!
 * \ingroup Geometry
 * \brief A class representing the intersection entities of two distributed geometric entity sets
 *
 * The intersections are distributed by domain ownership: each process holds exactly
 * those intersections whose domain (first) entity it owns. The union over all processes
 * equals the result of a sequential intersection of the global entity sets. The target
 * (second) entity of an intersection may be owned by a remote process; it is
 * identified by its global (rank, local index) identifier rather than by a local entity.
 *
 * \note This currently returns the raw (triangulated) intersections without the
 *       mixed-dimensional neighbor merging done by Dumux::IntersectionEntitySet.
 */
template<class DomainEntitySet, class TargetEntitySet>
class DistributedIntersectionEntitySet
{
    using DomainTree = DistributedBoundingBoxTree<DomainEntitySet>;
    using TargetTree = DistributedBoundingBoxTree<TargetEntitySet>;
    using Info = DistributedIntersectionInfo<DomainEntitySet, TargetEntitySet>;

public:
    using ctype = typename Info::ctype;
    static constexpr int dimensionworld = Info::dimensionworld;
    using GlobalPosition = typename Info::GlobalPosition;

private:
    static constexpr int dimDomain = DomainEntitySet::Entity::Geometry::mydimension;
    static constexpr int dimTarget = TargetEntitySet::Entity::Geometry::mydimension;
    static constexpr int dimIs = std::min(dimDomain, dimTarget);

public:
    //! the geometry of an intersection (always a simplex)
    using Geometry = Dune::AffineGeometry<ctype, dimIs, dimensionworld>;

    /*!
     * \brief A view on a single intersection entity
     */
    class IntersectionEntity
    {
    public:
        IntersectionEntity(const DistributedIntersectionEntitySet& set, const Info& info)
        : set_(set), info_(info) {}

        //! the intersection geometry
        Geometry geometry() const
        { return Geometry(Dune::GeometryTypes::simplex(dimIs), info_.corners); }

        //! the rank owning the domain (first) entity (always this process)
        int domainRank() const { return info_.domainRank; }
        //! the local index of the domain entity on its owner rank
        std::size_t domainIndex() const { return info_.domainIndex; }
        //! the rank owning the target (second) entity (possibly remote)
        int targetRank() const { return info_.targetRank; }
        //! the local index of the target entity on its owner rank
        std::size_t targetIndex() const { return info_.targetIndex; }

        //! the domain entity (always local to this process)
        typename DomainEntitySet::Entity domainEntity() const
        { return set_.domainTree_->entitySet().entity(info_.domainIndex); }

        //! whether the target entity is owned by this process
        bool targetIsLocal() const
        { return info_.targetRank == set_.targetTree_->rank(); }

        //! the target entity (only valid if targetIsLocal())
        typename TargetEntitySet::Entity targetEntity() const
        { return set_.targetTree_->entitySet().entity(info_.targetIndex); }

    private:
        const DistributedIntersectionEntitySet& set_;
        const Info& info_;
    };

    /*!
     * \brief An iterator over the intersection entities
     * \note dereferencing yields a lightweight IntersectionEntity view by value
     */
    class IntersectionIterator
    : public Dune::ForwardIteratorFacade<IntersectionIterator, const IntersectionEntity, IntersectionEntity>
    {
        using InfoIterator = typename std::vector<Info>::const_iterator;
    public:
        IntersectionIterator(const DistributedIntersectionEntitySet& set, InfoIterator it)
        : set_(&set), it_(it) {}

        IntersectionEntity dereference() const
        { return IntersectionEntity(*set_, *it_); }

        bool equals(const IntersectionIterator& other) const
        { return it_ == other.it_; }

        void increment() { ++it_; }

    private:
        const DistributedIntersectionEntitySet* set_;
        InfoIterator it_;
    };

    using EntityIterator = IntersectionIterator;

    //! Default constructor
    DistributedIntersectionEntitySet() = default;

    /*!
     * \brief Build the intersections from two distributed entity sets
     * \note collective operation over the processes of the entity sets
     */
    void build(std::shared_ptr<const DomainEntitySet> domainSet,
               std::shared_ptr<const TargetEntitySet> targetSet)
    {
        build(std::make_shared<DomainTree>(domainSet), std::make_shared<TargetTree>(targetSet));
    }

    /*!
     * \brief Build the intersections from two distributed bounding box trees
     * \note collective operation over the processes of the trees
     */
    void build(std::shared_ptr<const DomainTree> domainTree,
               std::shared_ptr<const TargetTree> targetTree)
    {
        domainTree_ = domainTree;
        targetTree_ = targetTree;

        Dune::Timer timer;
        intersections_ = intersectingEntities(*domainTree_, *targetTree_);

        const auto& comm = domainTree_->comm();
        const std::size_t globalSize = comm.sum(intersections_.size());
        if (comm.rank() == 0)
            std::cout << "Computed " << globalSize << " distributed intersection entities in "
                      << timer.elapsed() << " seconds." << std::endl;
    }

    //! the number of (process-local) intersections
    std::size_t size() const
    { return intersections_.size(); }

    //! begin iterator over the (process-local) intersections
    IntersectionIterator ibegin() const
    { return IntersectionIterator(*this, intersections_.begin()); }

    //! end iterator over the (process-local) intersections
    IntersectionIterator iend() const
    { return IntersectionIterator(*this, intersections_.end()); }

    /*!
     * \brief Range generator to iterate with range-based for loops over all intersections
     *        as follows: for (const auto& is : intersections(glue)) { ... }
     */
    friend Dune::IteratorRange<EntityIterator> intersections(const DistributedIntersectionEntitySet& set)
    { return {set.ibegin(), set.iend()}; }

private:
    std::vector<Info> intersections_;
    std::shared_ptr<const DomainTree> domainTree_;
    std::shared_ptr<const TargetTree> targetTree_;
};

} // end namespace Dumux

#endif
