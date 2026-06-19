// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Geometry
 * \brief Algorithms that find which geometric entities intersect across processes
 *        using a distributed (MPI-parallel) bounding box tree.
 */
#ifndef DUMUX_GEOMETRY_DISTRIBUTED_INTERSECTING_ENTITIES_HH
#define DUMUX_GEOMETRY_DISTRIBUTED_INTERSECTING_ENTITIES_HH

#include <vector>
#include <array>
#include <utility>

#include <dune/common/fvector.hh>
#include <dune/common/promotiontraits.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/partitionset.hh>

#include <dumux/parallel/packing.hh>
#include <dumux/parallel/nonblockingsparseexchange.hh>
#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/geometry/geometricentityset.hh>
#include <dumux/geometry/intersectingentities.hh>
#include <dumux/geometry/distributedboundingboxtree.hh>

namespace Dumux::Detail {

/*!
 * \brief Whether an entity belongs to one of the given parallel partitions
 * \note Entities without a partition concept (e.g. plain geometry wrappers) are
 *       always considered to belong: the only process that holds such an entity
 *       owns it.
 */
template<class Entity, class PartitionSet>
bool inPartition(const Entity& e, PartitionSet partitions)
{
    if constexpr (requires { e.partitionType(); })
        return partitions.contains(e.partitionType());
    else
        return true;
}

/*!
 * \brief Whether an entity is uniquely owned by the calling process
 * \note This is used for de-duplication: each entity of a distributed set is
 *       interior (and thus uniquely owned) on exactly one process. It must not
 *       be relaxed to include overlap/ghost entities, or the same entity would
 *       be contributed by more than one process.
 */
template<class Entity>
bool isUniquelyOwned(const Entity& e)
{ return inPartition(e, Dune::Partitions::interior); }

} // end namespace Dumux::Detail

namespace Dumux {

/*!
 * \ingroup Geometry
 * \brief Compute all intersections between entities and a point on a distributed tree
 * \param point the query point (typically identical on all processes)
 * \param tree the distributed bounding box tree
 * \param isCartesianGrid skip the primitive test for structured cube grids
 * \param onlyOwned only report entities owned (interior) by this process
 * \return the process-local intersections as (rank, local entity index) pairs
 *
 * \note This computes the intersections owned by the calling process only. Called
 *       collectively with the same point on all processes and merged, the union of
 *       results equals the result of a sequential query on the global entity set.
 */
template<class EntitySet, class ctype, int dimworld>
inline std::vector<std::pair<int, std::size_t>>
intersectingEntities(const Dune::FieldVector<ctype, dimworld>& point,
                     const DistributedBoundingBoxTree<EntitySet>& tree,
                     bool isCartesianGrid = false,
                     bool onlyOwned = true)
{
    std::vector<std::pair<int, std::size_t>> result;
    if (!tree.hasLocalEntities())
        return result;

    const auto localEntities = intersectingEntities(point, tree.localTree(), isCartesianGrid);
    const int rank = tree.rank();
    result.reserve(localEntities.size());
    for (const auto localIdx : localEntities)
    {
        const auto& set = tree.entitySet();
        if (onlyOwned && !Detail::isUniquelyOwned(set.entity(localIdx)))
            continue;
        result.emplace_back(rank, localIdx);
    }
    return result;
}

/*!
 * \ingroup Geometry
 * \brief An intersection of two entities of distributed entity sets
 * \note The intersecting entities are identified globally by a (rank, local index) pair
 */
template<class DomainEntitySet, class TargetEntitySet>
struct DistributedIntersectionInfo
{
    using ctype = typename Dune::PromotionTraits<typename DomainEntitySet::ctype, typename TargetEntitySet::ctype>::PromotedType;
    static constexpr int dimensionworld = DomainEntitySet::dimensionworld;
    using GlobalPosition = Dune::FieldVector<ctype, dimensionworld>;

    int domainRank; //!< the rank owning the domain (first) entity
    std::size_t domainIndex; //!< the local index of the domain entity on domainRank
    int targetRank; //!< the rank owning the target (second) entity
    std::size_t targetIndex; //!< the local index of the target entity on targetRank
    std::vector<GlobalPosition> corners; //!< the corners of the (simplex) intersection geometry
};

/*!
 * \ingroup Geometry
 * \brief Compute all intersections between two distributed bounding box trees
 * \param treeA the distributed bounding box tree of the domain entity set
 * \param treeB the distributed bounding box tree of the target entity set
 * \param domainPartitions the parallel partitions of the domain (treeA) entity set
 *        that this process should report intersections for (default: interior only)
 * \return the intersections whose domain entity is in domainPartitions on this process
 *
 * \note This is a collective operation. The result is distributed by domain partition:
 *       each process computes exactly the intersections whose domain (first) entity lies
 *       in domainPartitions on it. The target (second) entity is always identified by its
 *       unique (interior) owner so that no target entity is contributed more than once.
 *
 * \note With the default domainPartitions == Dune::Partitions::interior every domain entity
 *       is interior on exactly one process, so the union of all processes' results equals the
 *       result of a sequential query on the global entity sets (no duplicates). Pass e.g.
 *       Dune::Partitions::interior + Dune::Partitions::overlap to additionally report, on each
 *       process, the intersections of the domain entities it overlaps. This is what an
 *       overlapping (Schwarz) assembly needs — every process then holds the coupling
 *       intersections of all domain entities it assembles locally — at the price of overlap
 *       intersections being reported on more than one process (intended redundancy; not
 *       suitable for a global reduction over the result).
 *
 * Algorithm: each process broadcasts those of its uniquely-owned (interior) target entities
 * that overlap the partition bounding box of some other process. Every process then runs the
 * sequential tree-tree query of its local domain tree against (a) its own local target
 * tree and (b) a tree built over the imported target entities, keeping only those
 * intersections whose domain entity lies in domainPartitions and whose target entity it
 * uniquely owns (local pass) or imported from its owner (imported pass).
 */
template<class DomainEntitySet, class TargetEntitySet,
         class DomainPartitionSet = Dune::Partitions::Interior>
std::vector<DistributedIntersectionInfo<DomainEntitySet, TargetEntitySet>>
intersectingEntities(const DistributedBoundingBoxTree<DomainEntitySet>& treeA,
                     const DistributedBoundingBoxTree<TargetEntitySet>& treeB,
                     DomainPartitionSet domainPartitions = {})
{
    static_assert(int(DomainEntitySet::dimensionworld) == int(TargetEntitySet::dimensionworld),
        "Can only intersect distributed bounding box trees of the same world dimension");

    using Info = DistributedIntersectionInfo<DomainEntitySet, TargetEntitySet>;
    static constexpr int dimworld = DomainEntitySet::dimensionworld;
    static constexpr int dimTarget = TargetEntitySet::Entity::Geometry::mydimension;
    using TargetCtype = typename TargetEntitySet::ctype;
    using ImportedGeometry = Dune::MultiLinearGeometry<TargetCtype, dimTarget, dimworld>;
    using ImportedEntitySet = GeometriesEntitySet<ImportedGeometry>;
    using ImportedGlobalPosition = Dune::FieldVector<TargetCtype, dimworld>;

    const auto& comm = treeA.comm();
    const int myRank = comm.rank();
    const int numProc = comm.size();

    // 1. Pack each owned target entity into a buffer for every remote domain partition
    //    it overlaps. A single logarithmic tree-tree query of the local target tree
    //    against the replicated process tree of the domain determines, with subtree
    //    pruning, which target entities overlap which remote domain partitions.
    std::vector<std::vector<char>> sendBufs(numProc);
    if (treeB.hasLocalEntities())
    {
        const auto& setB = treeB.entitySet();

        // collect, per local target entity, the remote domain ranks it overlaps
        std::vector<std::vector<int>> destinations(setB.size());
        for (const auto& [targetIdx, processLeaf] : intersectingBoxes(treeB.localTree(), treeA.processTree()))
            if (const int p = treeA.processForLeaf(processLeaf); p != myRank)
                destinations[targetIdx].push_back(p);

        // pack each owned entity once and append it to each of its destinations' buffers
        std::vector<char> entityBuf;
        for (const auto& entity : setB)
        {
            const auto entityIdx = setB.index(entity);
            if (!Detail::isUniquelyOwned(entity) || destinations[entityIdx].empty())
                continue;

            entityBuf.clear();
            const auto geometry = entity.geometry();
            Detail::packValue(entityBuf, myRank);
            Detail::packValue(entityBuf, static_cast<std::size_t>(entityIdx));
            Detail::packValue(entityBuf, static_cast<unsigned int>(geometry.type().id()));
            const int numCorners = geometry.corners();
            Detail::packValue(entityBuf, numCorners);
            for (int c = 0; c < numCorners; ++c)
            {
                const auto corner = geometry.corner(c);
                for (int d = 0; d < dimworld; ++d)
                    Detail::packValue(entityBuf, static_cast<TargetCtype>(corner[d]));
            }

            for (const int p : destinations[entityIdx])
                sendBufs[p].insert(sendBufs[p].end(), entityBuf.begin(), entityBuf.end());
        }
    }

    // 2. Exchange the boundary target entities with their destination processes only.
    //    Each rank knows whom it sends to (sendBufs) but not who sends to it; a
    //    nonblocking-consensus (NBX) sparse data exchange discovers the senders and
    //    sizes on the fly, avoiding any dense O(numProc^2) count collective.
    std::vector<char> recvBuf;
    if (numProc > 1)
        recvBuf = Detail::exchangeSparse(comm, sendBufs);

    // 3. Deserialize imported target entities that overlap this process' domain box
    std::vector<ImportedGeometry> importedGeometries;
    std::vector<std::pair<int, std::size_t>> importedOrigin; // (rank, local index)
    if (treeA.hasLocalEntities() && !recvBuf.empty())
    {
        const auto* const myDomainBox = treeA.localBoundingBox().data();
        const char* cursor = recvBuf.data();
        const char* const end = recvBuf.data() + recvBuf.size();
        while (cursor < end)
        {
            const int originRank = Detail::unpackValue<int>(cursor);
            const auto originIndex = Detail::unpackValue<std::size_t>(cursor);
            const auto topologyId = Detail::unpackValue<unsigned int>(cursor);
            const int numCorners = Detail::unpackValue<int>(cursor);

            std::vector<ImportedGlobalPosition> corners(numCorners);
            for (int c = 0; c < numCorners; ++c)
                for (int d = 0; d < dimworld; ++d)
                    corners[c][d] = Detail::unpackValue<TargetCtype>(cursor);

            // skip our own entities (handled by the local-local pass below)
            if (originRank == myRank)
                continue;

            ImportedGeometry geometry(Dune::GeometryType(topologyId, dimTarget), corners);

            // skip entities that don't overlap our domain partition
            std::array<TargetCtype, 2*dimworld> box{};
            computeGeometryBoundingBox<dimworld>(box.data(), geometry);
            if (!intersectsBoundingBoxBoundingBox<dimworld>(box.data(), myDomainBox))
                continue;

            importedGeometries.push_back(std::move(geometry));
            importedOrigin.emplace_back(originRank, originIndex);
        }
    }

    // 4. Compute the local intersections, keeping only those we own (interior domain entity)
    std::vector<Info> result;
    if (!treeA.hasLocalEntities())
        return result;

    const auto& localTreeA = treeA.localTree();
    const auto& setA = treeA.entitySet();

    // (a) against our own local target entities
    if (treeB.hasLocalEntities())
    {
        const auto& setB = treeB.entitySet();
        const auto raw = intersectingEntities(localTreeA, treeB.localTree());
        for (const auto& is : raw)
        {
            if (!Detail::inPartition(setA.entity(is.first()), domainPartitions))
                continue;
            if (!Detail::isUniquelyOwned(setB.entity(is.second())))
                continue;
            result.push_back(Info{
                myRank, is.first(), myRank, is.second(),
                {is.corners().begin(), is.corners().end()}
            });
        }
    }

    // (b) against the imported (remote) target entities
    if (!importedGeometries.empty())
    {
        const auto importedSet = std::make_shared<const ImportedEntitySet>(std::move(importedGeometries));
        const BoundingBoxTree<ImportedEntitySet> importedTree(importedSet);
        const auto raw = intersectingEntities(localTreeA, importedTree);
        for (const auto& is : raw)
        {
            if (!Detail::inPartition(setA.entity(is.first()), domainPartitions))
                continue;
            const auto [originRank, originIndex] = importedOrigin[is.second()];
            result.push_back(Info{
                myRank, is.first(), originRank, originIndex,
                {is.corners().begin(), is.corners().end()}
            });
        }
    }

    return result;
}

} // end namespace Dumux

#endif
