// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedCoupling
 * \brief Helpers to set up the point sources of an embedded coupling on distributed grids
 */
#ifndef DUMUX_MULTIDOMAIN_EMBEDDED_DISTRIBUTED_POINT_SOURCES_HH
#define DUMUX_MULTIDOMAIN_EMBEDDED_DISTRIBUTED_POINT_SOURCES_HH

#include <algorithm>
#include <array>
#include <cstddef>
#include <span>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/geometry/intersectingentities.hh>
#include <dumux/geometry/distributedboundingboxtree.hh>
#include <dumux/geometry/distributedintersectingentities.hh>
#include <dumux/parallel/packing.hh>
#include <dumux/parallel/nonblockingsparseexchange.hh>

namespace Dumux::EmbeddedCoupling::Detail {

template<class T>
struct IsFieldVector : std::false_type {};

template<class T, int n>
struct IsFieldVector<Dune::FieldVector<T, n>> : std::true_type {};

template<class T>
struct IsPair : std::false_type {};

template<class T, class U>
struct IsPair<std::pair<T, U>> : std::true_type {};

//! Appends values, pairs, field vectors and vectors of trivially copyable types to a byte buffer
class MessageWriter
{
public:
    explicit MessageWriter(std::vector<char>& buffer) : buffer_(buffer) {}

    template<class T>
    void put(const T& value)
    {
        if constexpr (IsFieldVector<T>::value)
        {
            for (const auto& component : value)
                put(component);
        }
        else if constexpr (IsPair<T>::value)
        {
            put(value.first);
            put(value.second);
        }
        else
            Dumux::Detail::packValue(buffer_, value);
    }

    template<class T>
    void put(const std::vector<T>& values)
    {
        put(values.size());
        for (const auto& v : values)
            put(v);
    }

private:
    std::vector<char>& buffer_;
};

//! Reads what a MessageWriter wrote, in the same order
class MessageReader
{
public:
    explicit MessageReader(const std::vector<char>& buffer)
    : cursor_(buffer.data()), end_(buffer.data() + buffer.size()) {}

    template<class T>
    T get()
    {
        if constexpr (IsFieldVector<T>::value)
        {
            T value;
            for (auto& component : value)
                component = get<typename T::value_type>();
            return value;
        }
        else if constexpr (IsPair<T>::value)
        {
            auto first = get<typename T::first_type>();
            auto second = get<typename T::second_type>();
            return T(std::move(first), std::move(second));
        }
        else
            return Dumux::Detail::unpackValue<T>(cursor_);
    }

    template<class T>
    std::vector<T> getVector()
    {
        std::vector<T> values(get<std::size_t>());
        for (auto& v : values)
            v = get<T>();
        return values;
    }

    bool done() const
    { return cursor_ >= end_; }

private:
    const char* cursor_;
    const char* end_;
};

/*!
 * \brief Send one message to each of a sparse set of processes
 * \param comm the communication
 * \param messages one byte buffer per rank, empty if there is no message for that rank
 * \return the messages addressed to this process with their sender ranks, in ascending order of the
 *         sender rank; a message of this process to itself is returned as well
 * \note Collective over comm.
 */
template<class Communication>
std::vector<std::pair<int, std::vector<char>>>
exchangeMessages(const Communication& comm, std::vector<std::vector<char>>& messages, int tag = 7382)
{
    const int rank = comm.rank();
    std::vector<std::pair<int, std::vector<char>>> received;
    if (!messages[rank].empty())
        received.emplace_back(rank, std::move(messages[rank]));

    std::vector<std::vector<char>> sendBuffers(comm.size());
    for (int p = 0; p < comm.size(); ++p)
    {
        if (p == rank || messages[p].empty())
            continue;
        Dumux::Detail::packValue(sendBuffers[p], rank);
        Dumux::Detail::packValue(sendBuffers[p], messages[p].size());
        sendBuffers[p].insert(sendBuffers[p].end(), messages[p].begin(), messages[p].end());
    }

    const auto buffer = Dumux::Detail::exchangeSparse(comm, sendBuffers, tag);
    const char* cursor = buffer.data();
    const char* const end = buffer.data() + buffer.size();
    while (cursor < end)
    {
        const auto sender = Dumux::Detail::unpackValue<int>(cursor);
        const auto size = Dumux::Detail::unpackValue<std::size_t>(cursor);
        received.emplace_back(sender, std::vector<char>(cursor, cursor + size));
        cursor += size;
    }

    std::sort(received.begin(), received.end(), [](const auto& a, const auto& b){ return a.first < b.first; });
    return received;
}

/*!
 * \brief The ranks whose partition bounding box intersects the given box, in ascending order
 * \param tree the distributed bounding box tree
 * \param box the box, minimum then maximum coordinates
 * \param ranks the result
 * \param stack a buffer for the traversal
 */
template<class EntitySet, class ctype>
void ranksIntersectingBox(const DistributedBoundingBoxTree<EntitySet>& tree, const ctype* box,
                          std::vector<int>& ranks, std::vector<std::size_t>& stack)
{
    static constexpr int dimworld = EntitySet::dimensionworld;

    ranks.clear();
    const auto& processTree = tree.processTree();
    stack.assign(1, processTree.numBoundingBoxes() - 1);
    while (!stack.empty())
    {
        const auto nodeIdx = stack.back();
        stack.pop_back();
        if (!intersectsBoundingBoxBoundingBox<dimworld>(processTree.getBoundingBoxCoordinates(nodeIdx), box))
            continue;

        const auto& node = processTree.getBoundingBoxNode(nodeIdx);
        if (processTree.isLeaf(node, nodeIdx))
            ranks.push_back(tree.processForLeaf(node.child1));
        else
        {
            stack.push_back(node.child0);
            stack.push_back(node.child1);
        }
    }

    std::sort(ranks.begin(), ranks.end());
    ranks.erase(std::unique(ranks.begin(), ranks.end()), ranks.end());
}

/*!
 * \brief The ranks whose partition bounding box intersects the given box
 * \param tree the distributed bounding box tree
 * \param box the box, minimum then maximum coordinates
 */
template<class EntitySet, class ctype, std::size_t n>
std::vector<int> ranksIntersectingBox(const DistributedBoundingBoxTree<EntitySet>& tree,
                                      const std::array<ctype, n>& box)
{
    static_assert(n == 2*EntitySet::dimensionworld);
    std::vector<int> ranks;
    std::vector<std::size_t> stack;
    ranksIntersectingBox(tree, box.data(), ranks, stack);
    return ranks;
}

//! An interior element of a distributed grid: its rank, its index on that rank and optionally its geometry
template<class GlobalPosition>
struct LocatedElement
{
    int rank;
    std::size_t index;
    unsigned int topologyId = 0;
    std::vector<GlobalPosition> corners = {};

    friend bool operator<(const LocatedElement& a, const LocatedElement& b)
    { return std::tie(a.rank, a.index) < std::tie(b.rank, b.index); }
};

//! The interior elements containing each of a sequence of points
template<class GlobalPosition>
struct LocatedPoints
{
    std::vector<std::size_t> offsets = {0};
    std::vector<LocatedElement<GlobalPosition>> elements;

    //! the elements containing point i, sorted by rank and index
    std::span<const LocatedElement<GlobalPosition>> operator[](std::size_t i) const
    { return {elements.data() + offsets[i], offsets[i+1] - offsets[i]}; }

    std::size_t size() const
    { return offsets.size() - 1; }
};

//! Keeps all elements a search of the grid returns for a point
struct KeepAllElements
{
    template<class GlobalPosition>
    void operator()(std::vector<std::size_t>&, const GlobalPosition&) const {}
};

/*!
 * \brief For each point, the interior elements of a distributed grid that contain it
 * \param tree the distributed bounding box tree of the grid's elements
 * \param points the points of this process
 * \param withGeometry whether to return the topology id and the corners of the elements
 * \param removeElements removes elements from the indices of the elements of this process that the search
 *        returns, called with the indices and the point
 * \return per point all interior elements containing it, sorted by rank and index
 * \note Collective over the communication of the tree. The union of the elements over all ranks is the
 *       set a sequential query of the whole grid returns, without the removed elements. Only points in the
 *       partition boxes of other processes are sent to these processes.
 */
template<class EntitySet, class ctype, int dimworld, class RemoveElements = KeepAllElements>
LocatedPoints<Dune::FieldVector<ctype, dimworld>>
interiorElementsContaining(const DistributedBoundingBoxTree<EntitySet>& tree,
                           const std::vector<Dune::FieldVector<ctype, dimworld>>& points,
                           bool withGeometry = false,
                           const RemoveElements& removeElements = {})
{
    using GlobalPosition = Dune::FieldVector<ctype, dimworld>;
    using Located = LocatedElement<GlobalPosition>;
    const auto& comm = tree.comm();
    const int rank = comm.rank();

    // the interior elements of this process containing a point
    std::vector<std::size_t> candidates;
    const auto forLocalElementsContaining = [&](const GlobalPosition& point, auto&& f)
    {
        if (!tree.hasLocalEntities())
            return;
        const auto& localTree = tree.localTree();
        candidates.clear();
        intersectingEntities(point, localTree, localTree.numBoundingBoxes() - 1, candidates, false);
        std::erase_if(candidates, [&](const auto eIdx)
        { return !Dumux::Detail::isUniquelyOwned(tree.entitySet().entity(eIdx)); });
        removeElements(candidates, point);
        for (const auto eIdx : candidates)
            f(eIdx);
    };

    const auto makeLocated = [&](int owner, std::size_t eIdx, auto&& geometry)
    {
        Located element{owner, eIdx};
        if (withGeometry)
        {
            element.topologyId = geometry().type().id();
            element.corners.resize(geometry().corners());
            for (int c = 0; c < geometry().corners(); ++c)
                element.corners[c] = geometry().corner(c);
        }
        return element;
    };

    // the points in the partition box of this process are located here, the others are sent to the
    // processes whose partition boxes contain them
    LocatedPoints<GlobalPosition> local;
    local.offsets.reserve(points.size() + 1);
    std::vector<std::vector<char>> queries(comm.size());
    {
        std::vector<int> ranks;
        std::vector<std::size_t> stack;
        std::array<ctype, 2*dimworld> box;
        for (std::size_t i = 0; i < points.size(); ++i)
        {
            for (int d = 0; d < dimworld; ++d)
                box[d] = box[d + dimworld] = points[i][d];
            ranksIntersectingBox(tree, box.data(), ranks, stack);
            for (const int p : ranks)
            {
                if (p == rank)
                    forLocalElementsContaining(points[i], [&](std::size_t eIdx){
                        local.elements.push_back(makeLocated(rank, eIdx, [&]{ return tree.entitySet().entity(eIdx).geometry(); }));
                    });
                else
                {
                    MessageWriter writer(queries[p]);
                    writer.put(i);
                    writer.put(points[i]);
                }
            }
            std::sort(local.elements.begin() + local.offsets.back(), local.elements.end());
            local.offsets.push_back(local.elements.size());
        }
    }

    // answer with the interior elements of this process containing the points
    std::vector<std::vector<char>> answers(comm.size());
    for (const auto& [sender, message] : exchangeMessages(comm, queries))
    {
        MessageWriter writer(answers[sender]);
        MessageReader reader(message);
        while (!reader.done())
        {
            const auto i = reader.get<std::size_t>();
            const auto point = reader.get<GlobalPosition>();
            forLocalElementsContaining(point, [&](std::size_t eIdx){
                writer.put(i);
                writer.put(eIdx);
                if (withGeometry)
                {
                    const auto geometry = tree.entitySet().entity(eIdx).geometry();
                    writer.put(geometry.type().id());
                    writer.put(geometry.corners());
                    for (int c = 0; c < geometry.corners(); ++c)
                        writer.put(geometry.corner(c));
                }
            });
        }
    }

    std::vector<std::pair<std::size_t, Located>> remote;
    for (const auto& [sender, message] : exchangeMessages(comm, answers))
    {
        MessageReader reader(message);
        while (!reader.done())
        {
            const auto i = reader.get<std::size_t>();
            auto& element = remote.emplace_back(i, Located{sender, reader.get<std::size_t>()}).second;
            if (withGeometry)
            {
                element.topologyId = reader.get<unsigned int>();
                element.corners.resize(reader.get<int>());
                for (auto& corner : element.corners)
                    corner = reader.get<GlobalPosition>();
            }
        }
    }

    if (remote.empty())
        return local;

    // merge the elements of other processes into the elements of this process, point by point
    std::sort(remote.begin(), remote.end(), [](const auto& a, const auto& b){ return a.first < b.first; });
    LocatedPoints<GlobalPosition> result;
    result.offsets.reserve(points.size() + 1);
    result.elements.reserve(local.elements.size() + remote.size());
    auto next = remote.begin();
    for (std::size_t i = 0; i < points.size(); ++i)
    {
        const auto begin = result.elements.size();
        for (std::size_t j = local.offsets[i]; j < local.offsets[i+1]; ++j)
            result.elements.push_back(std::move(local.elements[j]));
        for (; next != remote.end() && next->first == i; ++next)
            result.elements.push_back(std::move(next->second));
        std::sort(result.elements.begin() + begin, result.elements.end());
        result.offsets.push_back(result.elements.size());
    }

    return result;
}

} // end namespace Dumux::EmbeddedCoupling::Detail

#endif
