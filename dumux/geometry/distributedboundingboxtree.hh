// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Geometry
 * \brief A distributed (MPI-parallel) axis-aligned bounding box volume hierarchy
 *
 * Each process builds a sequential Dumux::BoundingBoxTree over its local
 * partition of a distributed Dumux::GeometricEntitySet. The bounding box of
 * each local tree's root is communicated to all processes (allgather) and a
 * second, replicated bounding box tree (the "process tree") is built over
 * these per-process boxes. This upper level of the hierarchy is identical on
 * all processes and lets each process determine, with a logarithmic query and
 * without any further communication, which remote partitions a query might
 * collide with. The leaves of the process tree are the partitions, i.e. the
 * process-local trees, so the whole structure is a tree of trees.
 */
#ifndef DUMUX_GEOMETRY_DISTRIBUTED_BOUNDINGBOXTREE_HH
#define DUMUX_GEOMETRY_DISTRIBUTED_BOUNDINGBOXTREE_HH

#include <cassert>
#include <vector>
#include <array>
#include <span>
#include <memory>
#include <algorithm>
#include <type_traits>
#include <utility>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>

#include <dumux/geometry/boundingboxgeometry.hh>
#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/geometry/geometricentityset.hh>

namespace Dumux {

/*!
 * \ingroup Geometry
 * \brief An MPI-parallel axis-aligned bounding box volume tree
 *
 * A composition of a process-local Dumux::BoundingBoxTree and a replicated
 * Dumux::BoundingBoxTree over the per-process bounding boxes (the "process
 * tree"). The latter forms the upper, cross-process level of the hierarchy and
 * is identical on all ranks. Query it (processTree()) to find candidate processes
 * for a point or box, mapping the resulting leaves to ranks via processForLeaf(),
 * before querying the (process-local) tree of those candidates.
 *
 * \tparam GeometricEntitySet has the same requirements as for Dumux::BoundingBoxTree
 *         and additionally has to export a communicator via a comm() member
 *         function (e.g. Dumux::GridViewGeometricEntitySet).
 */
template<class GeometricEntitySet>
  requires requires(const GeometricEntitySet& set) { { set.comm() }; }
class DistributedBoundingBoxTree
{
    static constexpr int dimworld = GeometricEntitySet::dimensionworld;
    using ctype = typename GeometricEntitySet::ctype;

    // internal types for the replicated process tree (the tree over the per-process boxes,
    // each leaf being a process partition box carrying its owning rank)
    using ProcessBoxGeometry = BoundingBoxGeometry<ctype, dimworld>;
    using ProcessEntitySet = ProcessGeometricEntitySet<ctype, dimworld>;
    using ProcessEntity = typename ProcessEntitySet::Entity;

public:
    //! the type of the process-local bounding box tree
    using LocalTree = BoundingBoxTree<GeometricEntitySet>;
    //! the type of entity set this tree was built with
    using EntitySet = GeometricEntitySet;
    //! the communicator type
    using Communication = std::decay_t<decltype(std::declval<const GeometricEntitySet>().comm())>;
    //! the type of the replicated process tree
    using ProcessTree = BoundingBoxTree<ProcessEntitySet>;

    //! Default constructor
    DistributedBoundingBoxTree() = default;

    //! Constructor with an entity set
    explicit DistributedBoundingBoxTree(std::shared_ptr<const GeometricEntitySet> set)
    { build(set); }

    //! Build the distributed bounding box tree for the given entity set
    void build(std::shared_ptr<const GeometricEntitySet> set)
    {
        entitySet_ = set;
        comm_ = set->comm();
        const int numProc = comm_.size();

        // build the process-local tree (guarding against empty partitions for
        // which the sequential tree is not well-defined)
        const bool hasLocalEntities = set->size() > 0;
        if (hasLocalEntities)
            localTree_ = std::make_shared<const LocalTree>(set);
        else
            localTree_.reset();

        // determine the bounding box of the local partition (the root box of
        // the local tree) or a sentinel for empty partitions
        std::array<ctype, 2*dimworld> localBox{};
        if (hasLocalEntities)
            std::copy_n(localBoundingBox().data(), 2*dimworld, localBox.begin());

        // communicate the per-process bounding boxes to all processes
        std::vector<ctype> processBoxes(2*dimworld*numProc, 0.0);
        comm_.allgather(localBox.data(), 2*dimworld, processBoxes.data());

        // communicate which processes actually hold entities (valid boxes)
        const int localValid = hasLocalEntities ? 1 : 0;
        std::vector<int> processValid(numProc, 0);
        comm_.allgather(&localValid, 1, processValid.data());

        // build the replicated process tree over the non-empty process boxes;
        // each entity carries the rank it belongs to
        std::vector<ProcessEntity> processEntities;
        processEntities.reserve(numProc);
        for (int p = 0; p < numProc; ++p)
        {
            if (!processValid[p])
                continue;

            const ctype* b = processBoxes.data() + std::size_t(p)*2*dimworld;
            Dune::FieldVector<ctype, dimworld> lower, upper;
            for (int d = 0; d < dimworld; ++d)
            {
                lower[d] = b[d];
                upper[d] = b[d + dimworld];
            }
            // the entity's consecutive index is its position in the set
            processEntities.emplace_back(ProcessBoxGeometry{lower, upper}, processEntities.size(), p);
        }

        // the global entity set is assumed non-empty (at least one process holds
        // entities), just like for the sequential bounding box tree
        assert(!processEntities.empty());
        processTree_ = std::make_shared<const ProcessTree>(
            std::make_shared<const ProcessEntitySet>(std::move(processEntities)));
    }

    //! the entity set this tree was built with
    const EntitySet& entitySet() const
    { return *entitySet_; }

    //! the communicator this tree was built with
    const Communication& comm() const
    { return comm_; }

    //! the rank of this process
    int rank() const
    { return comm_.rank(); }

    //! whether this process holds any local entities
    bool hasLocalEntities() const
    { return static_cast<bool>(localTree_); }

    //! the process-local bounding box tree (only valid if hasLocalEntities())
    const LocalTree& localTree() const
    {
        if (!localTree_)
            DUNE_THROW(Dune::InvalidStateException, "No local bounding box tree on a process without entities");
        return *localTree_;
    }

    //! the replicated process tree over the per-process boxes
    const ProcessTree& processTree() const
    { return *processTree_; }

    //! the rank owning the given leaf of the process tree
    int processForLeaf(std::size_t leafIdx) const
    { return processTree_->entitySet().rank(leafIdx); }

    //! the bounding box (min then max) of this process' partition (only valid if hasLocalEntities())
    std::span<const ctype, 2*dimworld> localBoundingBox() const
    {
        if (!localTree_)
            DUNE_THROW(Dune::InvalidStateException, "No local bounding box on a process without entities");
        const auto* root = localTree_->getBoundingBoxCoordinates(localTree_->numBoundingBoxes() - 1);
        return std::span<const ctype, 2*dimworld>{root, 2*dimworld};
    }

private:
    //! the entity set this tree was built with
    std::shared_ptr<const EntitySet> entitySet_;

    //! the process-local bounding box tree (nullptr on empty partitions)
    std::shared_ptr<const LocalTree> localTree_;

    //! the replicated tree over the per-process boxes (nullptr if no process holds entities)
    std::shared_ptr<const ProcessTree> processTree_;

    //! the communicator
    Communication comm_{};
};

} // end namespace Dumux

#endif
