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
 *
 * \brief This class creates and manages the foreign overlap given an
 *        initial list of border indices and a BCRS matrix.
 *
 * The foreign overlap are all (row) indices which overlap with the
 * some of the current process's local indices.
 */
#ifndef DUMUX_FOREIGN_OVERLAP_FROM_BCRS_MATRIX_HH
#define DUMUX_FOREIGN_OVERLAP_FROM_BCRS_MATRIX_HH

#include <list>
#include <set>
#include <map>
#include <tuple>

#if HAVE_MPI
#include <mpi.h>
#endif // HAVE_MPI

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/operators.hh>

#include "borderindex.hh"

namespace Dumux {

/*!
 * \brief This class creates and manages the foreign overlap given an
 *        initial list of border indices and a BCRS matrix.
 *
 * The foreign overlap are all (row) indices which overlap with the
 * some of the current process's local indices.
 */
template<class BCRSMatrix>
class ForeignOverlapFromBCRSMatrix
{
    ForeignOverlapFromBCRSMatrix(const ForeignOverlapFromBCRSMatrix &A)
    {}

public:
    typedef int ProcessRank;
    typedef int BorderDistance;
    typedef int Index;
    typedef Index LocalIndex;
    typedef std::pair<LocalIndex, ProcessRank> IndexRank;
    typedef std::tuple<LocalIndex, ProcessRank, BorderDistance> IndexRankDist;
    typedef std::tuple<Index, BorderDistance, int> IndexDistanceNpeers;
    typedef std::list<IndexRankDist> SeedList;

    typedef std::set<ProcessRank> PeerSet;
    typedef std::list<IndexDistanceNpeers> ForeignOverlapWithPeer;
    typedef std::map<ProcessRank, ForeignOverlapWithPeer> ForeignOverlapByRank;
    typedef std::vector<std::map<ProcessRank, BorderDistance> > ForeignOverlapByIndex;

    typedef std::list<BorderIndex> BorderList;

    /*!
     * \brief Constructs the foreign overlap given a BCRS matrix and
     *        an initial list of border indices.
     */
    ForeignOverlapFromBCRSMatrix(const BCRSMatrix &A,
                                 const BorderList &foreignBorderList,
                                 const BorderList &domesticBorderList,
                                 int overlapSize)
        : foreignBorderList_(foreignBorderList)
        , domesticBorderList_(domesticBorderList)
    {
        overlapSize_ = overlapSize;

        myRank_ = 0;
#if HAVE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &myRank_);
#endif

        numLocal_ = A.N();

        // calculate the border list. From this, create an initial
        // seed list of indices which are in the overlap.
        SeedList initialSeedList;
        borderListToSeedList_(initialSeedList, foreignBorderList);

        // find the set of processes which have an overlap with the
        // local processes. (i.e. the set of processes which we will
        // have to communicate to.)
        seedListToPeerSet_(initialSeedList);

        // calculate the foreign overlap for the local partition,
        // i.e. find the distance of each row from the seed set.
        foreignOverlapByIndex_.resize(A.N());
        extendForeignOverlap_(A, initialSeedList, overlapSize);

        updateMasterRanks_(foreignBorderList, domesticBorderList);

        // group foreign overlap by peer process rank
        groupForeignOverlapByRank_();
    }

    /*!
     * \brief Returns true iff a local index is a border index.
     */
    bool isBorder(int localIdx) const
    { return borderIndices_.count(localIdx) > 0; }

    /*!
     * \brief Return the rank of the master process for a local index.
     */
    int masterOf(int localIdx) const
    { return masterRank_[localIdx]; }

    /*!
     * \brief Return true if the current rank is the "master" of an
     *        index.
     *
     * We define the master process as the process with the lowest
     * rank where the index is local. For a process is always the
     * master of its interior indices, but for border indices it is
     * only master if the index is not shared with a process of lower
     * rank.
     */
    bool iAmMasterOf(int localIdx) const
    { return masterRank_[localIdx] == myRank_; }

    /*!
     * \brief Returns the list of indices which intersect the process
     *        border and are in the interior of the local process.
     */
    const BorderList &foreignBorderList() const
    { return foreignBorderList_; }

    /*!
     * \brief Returns the list of indices which intersect the process
     *        border and are in the interior of some remote process.
     */
    const BorderList &domesticBorderList() const
    { return domesticBorderList_; }

    /*!
     * \brief Return true if a given local index is a domestic index
     *        for a peer.
     */
    bool isDomesticIndexFor(ProcessRank peerRank, Index localIdx) const
    {
        if (!isLocal(localIdx))
            // our own remote indices do not count!
            return false;

        typedef std::map<ProcessRank, BorderDistance> BorderDistMap;
        const BorderDistMap &borderDist = foreignOverlapByIndex_[localIdx];
        BorderDistMap::const_iterator bdIt = borderDist.find(peerRank);
        if (bdIt == borderDist.end())
            return false; // this index is not seen by the peer

        // the index is seen by the peer
        return true;
    }

    /*!
     * \brief Return the list of (local indices, border distance,
     *        number of processes) triples which are in the overlap of
     *        a given peer rank.
     */
    const ForeignOverlapWithPeer &foreignOverlapWithPeer(int peerRank) const
    {
        assert(foreignOverlapByRank_.find(peerRank) != foreignOverlapByRank_.end());
        return foreignOverlapByRank_.find(peerRank)->second;
    }

    /*!
     * \brief Return the map of (peer rank, border distance) for a given local index.
     */
    const std::map<ProcessRank, BorderDistance> &foreignOverlapByIndex(int localIdx) const
    {
        assert(isLocal(localIdx));
        return foreignOverlapByIndex_[localIdx];
    }

    /*!
     * \brief Return the set of process ranks which share an overlap
     *        with the current process.
     */
    const PeerSet &peerSet() const
    { return peerSet_; }

    /*!
     * \brief Returns the number local indices
     */
    int numLocal() const
    { return numLocal_; }

    /*!
     * \brief Returns true iff a domestic index is local
     */
    bool isLocal(int domesticIdx) const
    {
        return domesticIdx < numLocal();
    }

    /*!
     * \brief Returns true iff a local index is shared with an other process.
     */
    bool isShared(int domesticIdx) const
    { return isLocal(domesticIdx) && isShared_[domesticIdx]; }

    /*!
     * \brief Return the number of peer ranks for which a given local
     *        index is visible.
     */
    int numPeers(int localIdx) const
    { return foreignOverlapByIndex_[localIdx].size(); }

    /*!
     * \brief Returns the size of the overlap region
     */
    int overlapSize() const
    { return overlapSize_; }

    /*!
     * \brief Print the foreign overlap for debugging purposes.
     */
    void print() const
    {
        ForeignOverlapByRank::const_iterator it = foreignOverlapByRank_.begin();
        ForeignOverlapByRank::const_iterator endIt = foreignOverlapByRank_.end();
        for (; it != endIt; ++it) {
            std::cout << "Overlap rows(distance) for rank " << it->first << ": ";

            ForeignOverlapWithPeer::const_iterator rowIt = it->second.begin();
            ForeignOverlapWithPeer::const_iterator rowEndIt = it->second.end();
            for (; rowIt != rowEndIt; ++rowIt) {
                std::cout << std::get<0>(*rowIt) << "(" << std::get<1>(*rowIt) << ") ";
            }
            std::cout << "\n";
        }
    }

protected:
    // Calculate the set of peer processes from the initial seed list.
    void seedListToPeerSet_(const SeedList &seedList)
    {
        SeedList::const_iterator it = seedList.begin();
        SeedList::const_iterator endIt = seedList.end();
        for (; it != endIt; ++it)
            peerSet_.insert(std::get<1>(*it));
    }

    // calculate the local border indices given the initial seed list
    void borderListToSeedList_(SeedList &initialSeedList, const BorderList &borderList)
    {
        BorderList::const_iterator it = borderList.begin();
        BorderList::const_iterator endIt = borderList.end();
        for (; it != endIt; ++it) {
            initialSeedList.push_back(IndexRankDist(it->localIdx,
                                                    it->peerRank,
                                                    it->borderDistance));
            borderIndices_.insert(it->localIdx);
        }
    }

    // given a list of border indices and provided that
    // borderListToSeedList_() was already called, calculate the
    // master process of each local index.
    void updateMasterRanks_(const BorderList &foreignBorderList,
                            const BorderList &domesticBorderList)
    {
        // determine the minimum rank for all indices
        masterRank_.resize(numLocal_, myRank_);
        for (unsigned int localIdx = 0; localIdx < masterRank_.size(); ++localIdx) {
            int masterRank = myRank_;
            if (isBorder(localIdx)) {
                // if the local index is a border index, loop over all ranks
                // for which this index is also a border index. the lowest
                // rank wins!
                typedef typename std::map<ProcessRank, BorderDistance>::const_iterator iterator;
                iterator it = foreignOverlapByIndex_[localIdx].begin();
                iterator endIt = foreignOverlapByIndex_[localIdx].end();
                for (; it != endIt; ++it) {
                    if (it->second == 0) {
                        // if the border distance is zero, the rank with the minimum
                        masterRank = std::min(masterRank, it->first);
                    }
                }
            }
            masterRank_[localIdx] = masterRank;
        }

        // overwrite the master rank of the non-shared border indices
        isShared_.resize(numLocal_, false);
        BorderList::const_iterator it = foreignBorderList.begin();
        BorderList::const_iterator endIt = foreignBorderList.end();
        for (; it != endIt; ++it) {
            if (!it->isShared)
                masterRank_[it->localIdx] = myRank_;
            else
                isShared_[it->localIdx] = true;
        }

        // overwrite the master rank of the non-shared border on the
        // domestic overlap
        it = domesticBorderList_.begin();
        endIt = domesticBorderList_.end();
        for (; it != endIt; ++it) {
            if (!it->isShared) {
                masterRank_[it->localIdx] = it->peerRank;
            }
        }
    }


    // extend the foreign overlaps by one level. this uses a greedy
    // algorithm.
    void extendForeignOverlap_(const BCRSMatrix &A,
                               SeedList &seedList,
                               int overlapSize)
    {
        // add all processes in the seed rows of the current overlap
        // level
        int minOverlapDistance = overlapSize*2;
        SeedList::const_iterator it = seedList.begin();
        SeedList::const_iterator endIt = seedList.end();
        for (; it != endIt; ++it) {
            int localIdx = std::get<0>(*it);
            int peerRank = std::get<1>(*it);
            int distance = std::get<2>(*it);
            if (foreignOverlapByIndex_[localIdx].count(peerRank) == 0) {
                foreignOverlapByIndex_[localIdx][peerRank] = distance;
            }
            else {
                foreignOverlapByIndex_[localIdx][peerRank] =
                    std::min(distance,
                             foreignOverlapByIndex_[localIdx][peerRank]);
            }

            minOverlapDistance = std::min(minOverlapDistance, distance);
        }

        // if we have reached the maximum overlap distance, we're
        // finished and break the recursion
        if (minOverlapDistance >= overlapSize)
            return;

        // find the seed list for the next overlap level using the
        // seed set for the current level
        SeedList nextSeedList;
        it = seedList.begin();
        for (; it != endIt; ++it) {
            int rowIdx = std::get<0>(*it);
            int peerRank = std::get<1>(*it);
            int borderDist = std::get<2>(*it);

            // find all column indies in the row. The indices of the
            // columns are the additional indices of the overlap which
            // we would like to add
            typedef typename BCRSMatrix::ConstColIterator ColIterator;
            ColIterator colIt = A[rowIdx].begin();
            ColIterator colEndIt = A[rowIdx].end();
            for (; colIt != colEndIt; ++colIt) {
                int newIdx = colIt.index();

                // if the process is already is in the overlap of the
                // column index, ignore this column index!
                if (foreignOverlapByIndex_[newIdx].count(peerRank) > 0)
                    continue;


                // check whether the new index is already in the overlap
                bool hasIndex = false;
                typename SeedList::iterator sIt = nextSeedList.begin();
                typename SeedList::iterator sEndIt = nextSeedList.end();
                for (; sIt != sEndIt; ++sIt) {
                    if (std::get<0>(*sIt) == newIdx &&
                        std::get<1>(*sIt) == peerRank)
                    {
                        hasIndex = true;
                        std::get<2>(*sIt) = std::min(std::get<2>(*sIt), borderDist + 1);
                        break;
                    }
                }
                if (hasIndex)
                    continue; // we already have this index

                // add the current processes to the seed list for the
                // next overlap level
                IndexRankDist newTuple(newIdx, peerRank, borderDist + 1);
                nextSeedList.push_back(newTuple);
            }
        }

        // clear the old seed list to save some memory
        seedList.clear();

        // Perform the same excercise for the next overlap distance
        extendForeignOverlap_(A,
                              nextSeedList,
                              overlapSize);
    }

    // assuming that the foreign overlap has been created for each
    // local index, this method groups the foreign overlap by peer
    // process rank
    void groupForeignOverlapByRank_()
    {
        // loop over all indices which are in the overlap of some
        // process
        int nIndices = foreignOverlapByIndex_.size();
        for (int i = 0; i < nIndices; ++i)
        {
            // loop over the list of processes for the current index
            std::map<ProcessRank, BorderDistance>::const_iterator it =
                foreignOverlapByIndex_[i].begin();
            std::map<ProcessRank, BorderDistance>::const_iterator endIt =
                foreignOverlapByIndex_[i].end();
            int nRanks = foreignOverlapByIndex_[i].size();
            for (; it != endIt; ++it)  {
                IndexDistanceNpeers tmp(i, it->second, nRanks);
                foreignOverlapByRank_[it->first].push_back(tmp);
            }
        }
    }

    // set of processes with which we have to communicate
    PeerSet peerSet_;

    // the list of indices on the border
    const BorderList &foreignBorderList_;
    const BorderList &domesticBorderList_;

    // an array which contains the rank of the master process for each
    // index
    std::vector<ProcessRank> masterRank_;

    // an array which stores whether an index is also in the interior
    // of some other process
    std::vector<bool> isShared_;

    // set of all local indices which are on the border of some remote
    // process
    std::set<LocalIndex> borderIndices_;

    // stores the set of process ranks which are in the overlap for a
    // given row index "owned" by the current rank. The second value
    // store the distance from the nearest process border.
    ForeignOverlapByIndex foreignOverlapByIndex_;

    // stores a list of foreign overlap indices for each rank
    ForeignOverlapByRank foreignOverlapByRank_;

    // extend of the overlap region
    int overlapSize_;

    // number of local indices
    int numLocal_;

    // the MPI rank of the local process
    int myRank_;
};

} // namespace Dumux

#endif
