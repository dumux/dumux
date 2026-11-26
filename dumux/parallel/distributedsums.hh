// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Parallel
 * \brief Sums whose terms are computed on several processes
 */
#ifndef DUMUX_PARALLEL_DISTRIBUTED_SUMS_HH
#define DUMUX_PARALLEL_DISTRIBUTED_SUMS_HH

#include <algorithm>
#include <cstddef>
#include <limits>
#include <map>
#include <numeric>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>

#if HAVE_MPI
#include <mpi.h>
#endif

namespace Dumux {

/*!
 * \ingroup Parallel
 * \brief Sums whose terms are computed on several processes
 *
 * Every sum has a global key, the rank that created it and an index on that rank, and a fixed number of
 * scalar components. A process holding a term of a sum sends it to the processes that need the sum,
 * which add up the terms they receive. The pattern of who sends which term to whom is fixed by
 * finalize(); every exchange() is then one round of point-to-point messages between neighboring
 * processes, each holding the terms in the order of the keys.
 *
 * \note A process does not send to itself: its own term is not part of received().
 */
template<class Scalar>
class DistributedSums
{
public:
    //! the rank that created a sum and the index of the sum on that rank
    using Key = std::pair<int, std::size_t>;

    static constexpr std::size_t invalidIndex = std::numeric_limits<std::size_t>::max();

    DistributedSums() = default;

    template<class Communication>
    explicit DistributedSums(const Communication& comm)
    : rank_(comm.rank())
    , size_(comm.size())
    {
#if HAVE_MPI
        if constexpr (std::is_convertible_v<Communication, MPI_Comm>)
            comm_ = comm;
#endif
    }

    //! The local index of the sum with the given key, added with numComponents components if it is new
    std::size_t insert(const Key& key, std::size_t numComponents)
    {
        if (finalized_)
            DUNE_THROW(Dune::InvalidStateException, "Cannot add sums after finalize()");

        const auto [it, inserted] = keyToIndex_.try_emplace(key, keys_.size());
        if (inserted)
        {
            keys_.push_back(key);
            offsets_.push_back(offsets_.back() + numComponents);
        }
        else if (this->numComponents(it->second) != numComponents)
            DUNE_THROW(Dune::InvalidStateException, "A sum was added again with a different number of components");

        return it->second;
    }

    //! The local index of the sum with the given key, invalidIndex if this process does not know it
    std::size_t index(const Key& key) const
    {
        const auto it = keyToIndex_.find(key);
        return it == keyToIndex_.end() ? invalidIndex : it->second;
    }

    //! The key of the sum with local index i
    const Key& key(std::size_t i) const
    { return keys_[i]; }

    //! The number of sums known to this process
    std::size_t size() const
    { return keys_.size(); }

    //! The number of components of the sum with local index i
    std::size_t numComponents(std::size_t i) const
    { return offsets_[i+1] - offsets_[i]; }

    //! This process holds a term of sum i and sends it to the given rank
    void addDestination(std::size_t i, int rank)
    {
        if (rank != rank_)
            destinations_.emplace_back(rank, i);
    }

    //! This process needs sum i and receives the term held by the given rank
    void addSource(std::size_t i, int rank)
    {
        if (rank != rank_)
            sources_.emplace_back(rank, i);
    }

    /*!
     * \brief Fix the communication pattern
     * \note Not collective. The pattern has to be consistent: if rank p adds rank q as a destination of a
     *       sum, rank q adds rank p as a source of the same sum.
     */
    void finalize()
    {
        makePeers_(destinations_, sendPeers_);
        makePeers_(sources_, receivePeers_);
        receives_.assign(size(), false);
        for (const auto& peer : receivePeers_)
            for (const auto i : peer.sums)
                receives_[i] = true;
        received_.assign(offsets_.back(), 0.0);
        finalized_ = true;
    }

    /*!
     * \brief Exchange the terms and add up the received ones
     * \param term the function term(i, values) writing this process' term of sum i into the span values;
     *        called for every sum this process sends
     * \note Collective over the processes that are peers of this process.
     */
    template<class Term>
    void exchange(Term&& term)
    { exchange(std::forward<Term>(term), received_); }

    /*!
     * \brief Exchange the terms and add up the received ones into the given values instead of received()
     * \param term the function term(i, values) writing this process' term of sum i into the span values;
     *        called for every sum this process sends
     * \param values the received sums, read with received(i, values)
     * \note Collective over the processes that are peers of this process.
     */
    template<class Term>
    void exchange(Term&& term, std::vector<Scalar>& values)
    {
        if (!finalized_)
            DUNE_THROW(Dune::InvalidStateException, "Call finalize() before exchange()");

        values.assign(offsets_.back(), 0.0);
        if (sendPeers_.empty() && receivePeers_.empty())
            return;

#if HAVE_MPI
        std::vector<MPI_Request> requests;
        requests.reserve(sendPeers_.size() + receivePeers_.size());
        for (auto& peer : receivePeers_)
        {
            peer.buffer.resize(peer.numScalars);
            requests.emplace_back();
            MPI_Irecv(peer.buffer.data(), static_cast<int>(peer.numScalars*sizeof(Scalar)), MPI_BYTE,
                      peer.rank, tag_, comm_, &requests.back());
        }

        for (auto& peer : sendPeers_)
        {
            peer.buffer.resize(peer.numScalars);
            std::size_t offset = 0;
            for (const auto i : peer.sums)
            {
                const auto n = numComponents(i);
                std::fill_n(peer.buffer.begin() + offset, n, 0.0);
                term(i, std::span<Scalar>(peer.buffer.data() + offset, n));
                offset += n;
            }
            requests.emplace_back();
            MPI_Isend(peer.buffer.data(), static_cast<int>(peer.numScalars*sizeof(Scalar)), MPI_BYTE,
                      peer.rank, tag_, comm_, &requests.back());
        }

        std::vector<MPI_Status> statuses(requests.size());
        MPI_Waitall(static_cast<int>(requests.size()), requests.data(), statuses.data());

        for (std::size_t p = 0; p < receivePeers_.size(); ++p)
        {
            const auto& peer = receivePeers_[p];
            int count = 0;
            MPI_Get_count(&statuses[p], MPI_BYTE, &count);
            if (static_cast<std::size_t>(count) != peer.numScalars*sizeof(Scalar))
                DUNE_THROW(Dune::InvalidStateException, "Rank " << rank_ << " expected " << peer.numScalars
                           << " values from rank " << peer.rank << " but received " << count/sizeof(Scalar));

            std::size_t offset = 0;
            for (const auto i : peer.sums)
                for (std::size_t c = 0; c < numComponents(i); ++c)
                    values[offsets_[i] + c] += peer.buffer[offset++];
        }
#else
        DUNE_THROW(Dune::InvalidStateException, "Distributed sums with peer processes need MPI");
#endif
    }

    //! The sum of the terms of sum i received from other processes in the last exchange
    std::span<const Scalar> received(std::size_t i) const
    { return received(i, received_); }

    //! The sum of the terms of sum i in values filled by exchange(term, values)
    std::span<const Scalar> received(std::size_t i, const std::vector<Scalar>& values) const
    { return std::span<const Scalar>(values.data() + offsets_[i], numComponents(i)); }

    //! Whether this process receives terms of sum i from other processes
    bool receives(std::size_t i) const
    { return receives_[i]; }

    /*!
     * \brief Send lists of values along the pattern of the sums
     * \param list the function list(i) returning the std::vector<T> of this process for sum i; called for
     *        every sum this process sends
     * \return for every sum this process receives, the lists of all processes it receives from, in the order
     *         of their ranks; empty for the other sums
     * \note Collective over the processes that are peers of this process. T has to be trivially copyable.
     */
    template<class T, class List>
    std::vector<std::vector<T>> exchangeLists(List&& list) const
    {
        static_assert(std::is_trivially_copyable_v<T>, "The values of the lists are sent as bytes");
        if (!finalized_)
            DUNE_THROW(Dune::InvalidStateException, "Call finalize() before exchangeLists()");

        std::vector<std::vector<T>> lists(size());
        if (sendPeers_.empty() && receivePeers_.empty())
            return lists;

#if HAVE_MPI
        std::vector<std::vector<std::size_t>> sendCounts(sendPeers_.size());
        std::vector<std::vector<T>> sendValues(sendPeers_.size());
        for (std::size_t p = 0; p < sendPeers_.size(); ++p)
        {
            for (const auto i : sendPeers_[p].sums)
            {
                const std::vector<T> values = list(i);
                sendCounts[p].push_back(values.size());
                sendValues[p].insert(sendValues[p].end(), values.begin(), values.end());
            }
        }

        std::vector<std::vector<std::size_t>> receiveCounts(receivePeers_.size());
        for (std::size_t p = 0; p < receivePeers_.size(); ++p)
            receiveCounts[p].resize(receivePeers_[p].sums.size());
        sendReceive_(sendCounts, receiveCounts, listTag_);

        std::vector<std::vector<T>> receiveValues(receivePeers_.size());
        for (std::size_t p = 0; p < receivePeers_.size(); ++p)
            receiveValues[p].resize(std::accumulate(receiveCounts[p].begin(), receiveCounts[p].end(), std::size_t(0)));
        sendReceive_(sendValues, receiveValues, listTag_ + 1);

        for (std::size_t p = 0; p < receivePeers_.size(); ++p)
        {
            std::size_t offset = 0;
            const auto& sums = receivePeers_[p].sums;
            for (std::size_t k = 0; k < sums.size(); ++k)
            {
                auto& target = lists[sums[k]];
                target.insert(target.end(), receiveValues[p].begin() + offset, receiveValues[p].begin() + offset + receiveCounts[p][k]);
                offset += receiveCounts[p][k];
            }
        }
#else
        DUNE_THROW(Dune::InvalidStateException, "Distributed sums with peer processes need MPI");
#endif
        return lists;
    }

    //! The ranks this process sends its terms to, in ascending order
    std::vector<int> destinationRanks() const
    { return ranks_(sendPeers_); }

    //! The ranks this process receives terms from, in ascending order
    std::vector<int> sourceRanks() const
    { return ranks_(receivePeers_); }

private:
    struct Peer
    {
        int rank;
        std::vector<std::size_t> sums;
        std::size_t numScalars = 0;
        std::vector<Scalar> buffer;
    };

    void makePeers_(std::vector<std::pair<int, std::size_t>>& links, std::vector<Peer>& peers) const
    {
        std::sort(links.begin(), links.end(), [&](const auto& a, const auto& b)
        { return a.first != b.first ? a.first < b.first : keys_[a.second] < keys_[b.second]; });
        links.erase(std::unique(links.begin(), links.end()), links.end());

        if (!links.empty() && size_ == 1)
            DUNE_THROW(Dune::InvalidStateException, "Distributed sums with peer processes need a parallel communicator");

        peers.clear();
        for (const auto& [rank, i] : links)
        {
            if (peers.empty() || peers.back().rank != rank)
                peers.push_back(Peer{rank, {}, 0, {}});
            peers.back().sums.push_back(i);
            peers.back().numScalars += numComponents(i);
        }
    }

#if HAVE_MPI
    //! Send send[p] to the p-th send peer and receive receive[p], whose size is known, from the p-th receive peer
    template<class T>
    void sendReceive_(const std::vector<std::vector<T>>& send, std::vector<std::vector<T>>& receive, int tag) const
    {
        std::vector<MPI_Request> requests;
        requests.reserve(sendPeers_.size() + receivePeers_.size());
        for (std::size_t p = 0; p < receivePeers_.size(); ++p)
        {
            requests.emplace_back();
            MPI_Irecv(receive[p].data(), static_cast<int>(receive[p].size()*sizeof(T)), MPI_BYTE,
                      receivePeers_[p].rank, tag, comm_, &requests.back());
        }
        for (std::size_t p = 0; p < sendPeers_.size(); ++p)
        {
            requests.emplace_back();
            MPI_Isend(send[p].data(), static_cast<int>(send[p].size()*sizeof(T)), MPI_BYTE,
                      sendPeers_[p].rank, tag, comm_, &requests.back());
        }
        MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
    }
#endif

    static std::vector<int> ranks_(const std::vector<Peer>& peers)
    {
        std::vector<int> ranks;
        for (const auto& peer : peers)
            ranks.push_back(peer.rank);
        return ranks;
    }

    int rank_ = 0;
    int size_ = 1;
#if HAVE_MPI
    MPI_Comm comm_ = MPI_COMM_SELF;
#endif
    static constexpr int tag_ = 7381;
    static constexpr int listTag_ = 7382;

    std::vector<Key> keys_;
    std::vector<std::size_t> offsets_ = {0};
    std::map<Key, std::size_t> keyToIndex_;

    std::vector<std::pair<int, std::size_t>> destinations_, sources_;
    std::vector<Peer> sendPeers_, receivePeers_;
    std::vector<Scalar> received_;
    std::vector<bool> receives_;
    bool finalized_ = false;
};

} // end namespace Dumux

#endif
