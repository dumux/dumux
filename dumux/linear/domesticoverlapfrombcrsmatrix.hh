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
#ifndef DUMUX_DOMESTIC_OVERLAP_FROM_BCRS_MATRIX_HH
#define DUMUX_DOMESTIC_OVERLAP_FROM_BCRS_MATRIX_HH

#warning This file is deprecated and will be removed after Dumux 2.9

#include <list>
#include <set>
#include <map>
#include <memory>
#include <tuple>

#include <dumux/parallel/mpibuffer.hh>

#include "foreignoverlapfrombcrsmatrix.hh"
#include "globalindices.hh"

namespace Dumux {

/*!
 * \brief This class creates and manages the foreign overlap given an
 *        initial list of border indices and a BCRS matrix.
 *
 * The foreign overlap are all (row) indices which overlap with the
 * some of the current process's local indices.
 */
template<class BCRSMatrix>
class DomesticOverlapFromBCRSMatrix
{
    DomesticOverlapFromBCRSMatrix(const DomesticOverlapFromBCRSMatrix &A)
    {}

    typedef Dumux::ForeignOverlapFromBCRSMatrix<BCRSMatrix> ForeignOverlap;
    typedef Dumux::GlobalIndices<ForeignOverlap> GlobalIndices;

public:
    typedef int ProcessRank;
    typedef int BorderDistance;
    typedef int Index;
    typedef Index LocalIndex;
    typedef std::tuple<Index, BorderDistance, int> IndexDistanceNpeers;

    typedef std::set<ProcessRank> PeerSet;

    typedef std::list<IndexDistanceNpeers> ForeignOverlapWithPeer;

    typedef std::pair<LocalIndex, BorderDistance> IndexDistance;
    typedef std::set<LocalIndex> DomesticOverlapWithPeer;
    typedef std::map<ProcessRank, DomesticOverlapWithPeer> DomesticOverlapByRank;
    typedef std::vector<std::map<ProcessRank, BorderDistance> > DomesticOverlapByIndex;

    typedef std::list<BorderIndex> BorderList;

    /*!
     * \brief Constructs the foreign overlap given a BCRS matrix and
     *        an initial list of border indices.
     */
    DomesticOverlapFromBCRSMatrix(const BCRSMatrix &A,
                                  const BorderList &borderList,
                                  const BorderList &domesticBorderList,
                                  int overlapSize)
        : foreignOverlap_(A, borderList, domesticBorderList, overlapSize)
        , globalIndices_(foreignOverlap_)
    {
        myRank_ = 0;

#if HAVE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &myRank_);
#endif // HAVE_MPI

        buildDomesticOverlap_(A);
    }

    /*!
     * \brief Return the border list.
     *
     * The border list is the list of (local index, peer index, peer
     * rank) triples for all indices on a process border.
     */
    const BorderList& borderList() const
    { return foreignOverlap_.borderList(); }


    /*!
     * \brief Returns true iff a domestic index is a border index.
     */
    bool isBorder(int domesticIdx) const
    {
        return isLocal(domesticIdx) && foreignOverlap_.isBorder(domesticIdx);
    }

    /*!
     * \brief Returns true iff a domestic index is a front index.
     */
    bool isFront(int domesticIdx) const
    { return borderDistance_[domesticIdx] == foreignOverlap_.overlapSize(); }

    /*!
     * \brief Returns the number of processes which "see" a given
     *        index.
     */
    int numPeers(int domesticIdx) const
    { return domesticOverlapByIndex_[domesticIdx].size(); }

    /*!
     * \brief Returns whether a given domestic index is a front index
     *        for a given peer process.
     */
    int isFrontFor(int peerRank, int domesticIdx) const
    {
        typedef std::map<ProcessRank, BorderDistance>::const_iterator iterator;
        iterator it = domesticOverlapByIndex_[domesticIdx].find(peerRank);
        if (it == domesticOverlapByIndex_[domesticIdx].end())
            return false; // not seen by the process

        return it->second == foreignOverlap_.overlapSize();
    }

    /*!
     * \brief Return the number of processes which "see" a domestic
     *        index and for which the index is not on the front.
     */
    int numNonFrontProcesses(int domesticIdx) const
    {
        int result = 0;
        if (!isFront(domesticIdx))
            ++result;

        typedef std::map<ProcessRank, BorderDistance>::const_iterator iterator;
        iterator it = domesticOverlapByIndex_[domesticIdx].begin();
        iterator endIt = domesticOverlapByIndex_[domesticIdx].end();
        for (; it != endIt; ++it) {
            if (it->second < overlapSize())
                ++result;
        }

        return result;
    }

    /*!
     * \brief Returns the rank of the current process.
     */
    int myRank() const
    { return myRank_; }

    /*!
     * \brief Return the set of process ranks which share an overlap
     *        with the current process.
     */
    const PeerSet &peerSet() const
    { return peerSet_; }


    /*!
     * \brief Returns the foreign overlap
     */
    const ForeignOverlap &foreignOverlap() const
    { return foreignOverlap_; }

    /*!
     * \brief Returns the size of the overlap region
     */
    int overlapSize() const
    { return foreignOverlap_.overlapSize(); }

    /*!
     * \brief Returns the foreign overlap of a peer process
     */
    const ForeignOverlapWithPeer &foreignOverlapWithPeer(ProcessRank peerRank) const
    { return foreignOverlap_.foreignOverlapWithPeer(peerRank); }

    /*!
     * \brief Returns the domestic overlap with a peer process
     */
    const DomesticOverlapWithPeer &domesticOverlapWithPeer(ProcessRank peerRank) const
    {
        assert(domesticOverlapWithPeer_.find(peerRank) != domesticOverlapWithPeer_.end());
        return domesticOverlapWithPeer_.find(peerRank)->second;
    }

    /*!
     * \brief Returns true iff a given local index is a remote index for a given peer
     */
    bool isRemoteIndexFor(ProcessRank peerRank, Index localIdx) const
    { return foreignOverlap_.isRemoteIndexFor(peerRank, localIdx); }

    /*!
     * \brief Returns true iff a given local index is also a local index for a given peer
     */
    bool isLocalIndexFor(ProcessRank peerRank, Index domesticIdx) const
    { return foreignOverlap_.isLocalIndexFor(peerRank, domesticIdx); }

    /*!
     * \brief Returns true iff a given local index is a domestic index for a given peer
     */
    bool isDomesticIndexFor(ProcessRank peerRank, Index domesticIdx) const
    { return foreignOverlap_.isDomesticIndexFor(peerRank, domesticIdx); }

    /*!
     * \brief Returns the number local indices
     */
    int numLocal() const
    { return foreignOverlap_.numLocal(); }

    /*!
     * \brief Returns the number domestic indices.
     *
     * The domestic indices are defined as the process' local indices
     * plus its copies of indices in the overlap regions
     */
    int numDomestic() const
    { return globalIndices_.numDomestic(); }

    /*!
     * \brief Return true if a domestic index is local for the process
     *        (i.e. interior or border)
     */
    bool isLocal(int domesticIdx) const
    { return domesticIdx < numLocal(); }

    /*!
     * \brief Return true iff the current process is the master of a
     *        given domestic index.
     */
    bool iAmMasterOf(int domesticIdx) const
    {
        if (!isLocal(domesticIdx))
            return false;
        return foreignOverlap_.masterOf(domesticIdx) == myRank_;
    }

    /*!
     * \brief Return true iff a given index is shared by more than one process
     */
    bool isShared(int domesticIdx) const
    {
        if (!isLocal(domesticIdx))
            return false;
        return foreignOverlap_.isShared(domesticIdx);
    }

    /*!
     * \brief Return true iff a given rank is the master of a given
     *        domestic index.
     */
    bool isMasterOf(int peerRank, int domesticIdx) const
    {
        if (isLocal(domesticIdx)) {
            return foreignOverlap_.masterOf(domesticIdx) == peerRank;
        }

        // if the local index is a border index, loop over all ranks
        // for which this index is also a border index. the lowest
        // rank wins!
        typedef typename std::map<ProcessRank, BorderDistance>::const_iterator iterator;
        iterator it = domesticOverlapByIndex_[domesticIdx].begin();
        iterator endIt = domesticOverlapByIndex_[domesticIdx].end();
        LocalIndex masterIdx = std::numeric_limits<int>::max();
        for (; it != endIt; ++it) {
            masterIdx = std::min(masterIdx, it->first);
        }

        return masterIdx == peerRank;
    }

    /*!
     * \brief Print the foreign overlap for debugging purposes.
     */
    void print() const
    {
        globalIndices_.print();
    }

    /*!
     * \brief Returns a domestic index given a global one
     */
    Index globalToDomestic(Index globalIdx) const
    {
        return globalIndices_.globalToDomestic(globalIdx);
    }

    /*!
     * \brief Returns a global index given a domestic one
     */
    Index domesticToGlobal(Index domIdx) const
    {
        return globalIndices_.domesticToGlobal(domIdx);
    }

protected:
    void buildDomesticOverlap_(const BCRSMatrix &A)
    {
        // copy the set of peers from the foreign overlap
        peerSet_ = foreignOverlap_.peerSet();

        // resize the array which stores the number of peers for
        // each entry.
        domesticOverlapByIndex_.resize(numLocal());
        borderDistance_.resize(numLocal(), -1);

        // for all local indices copy the number of processes from the
        // foreign overlap
        for (int localIdx = 0; localIdx < numLocal(); ++localIdx) {
            const std::map<ProcessRank, BorderDistance> &idxOverlap =
                foreignOverlap_.foreignOverlapByIndex(localIdx);
            std::map<ProcessRank, BorderDistance>::const_iterator it = idxOverlap.begin();
            std::map<ProcessRank, BorderDistance>::const_iterator endIt = idxOverlap.end();
            for (; it != endIt; ++it) {
                int peerRank = it->first;
                int borderDist = it->second;
                domesticOverlapByIndex_[localIdx][peerRank] = borderDist;
                domesticOverlapWithPeer_[peerRank].insert(localIdx);
            }
        }

        PeerSet::const_iterator peerIt;
        PeerSet::const_iterator peerEndIt = foreignOverlap_.peerSet().end();

        // send the overlap indices to all peer processes
        peerIt = foreignOverlap_.peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            sendIndicesToPeer_(peerRank);
        }

        // receive our overlap from the processes to all peer processes
        peerIt = foreignOverlap_.peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            receiveIndicesFromPeer_(peerRank);
        }

        // receive our overlap from the processes to all peer processes
        peerIt = foreignOverlap_.peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            waitSendIndices_(peerRank);
        }
    }

    void sendIndicesToPeer_(int peerRank)
    {
#if HAVE_MPI
        const ForeignOverlapWithPeer &foreignOverlap
            = foreignOverlap_.foreignOverlapWithPeer(peerRank);

        // first, send a message containing the number of additional
        // indices stemming from the overlap (i.e. without the border
        // indices)
        int numIndices = foreignOverlap.size();
        numIndicesSendBuff_[peerRank] = std::make_shared<MpiBuffer<int> >(1);
        (*numIndicesSendBuff_[peerRank])[0] = numIndices;
        numIndicesSendBuff_[peerRank]->send(peerRank);

        // create MPI buffers
        indicesSendBuff_[peerRank] = std::make_shared<MpiBuffer<IndexDistanceNpeers> >(numIndices);

        // then send the additional indices themselfs
        ForeignOverlapWithPeer::const_iterator overlapIt = foreignOverlap.begin();
        ForeignOverlapWithPeer::const_iterator overlapEndIt = foreignOverlap.end();
        for (int i = 0; overlapIt != overlapEndIt; ++overlapIt, ++i) {
            int localIdx = std::get<0>(*overlapIt);
            int borderDistance = std::get<1>(*overlapIt);

            const std::map<ProcessRank, BorderDistance> &foreignIndexOverlap
                = foreignOverlap_.foreignOverlapByIndex(localIdx);

            int numPeers = foreignIndexOverlap.size();

            (*indicesSendBuff_[peerRank])[i] =
                IndexDistanceNpeers(globalIndices_.domesticToGlobal(localIdx),
                                    borderDistance,
                                    numPeers);

            // send all peer ranks which see the given index
            peersSendBuff_[peerRank].push_back(std::make_shared<MpiBuffer<int> >(2*numPeers));
            typename std::map<ProcessRank, BorderDistance>::const_iterator it = foreignIndexOverlap.begin();
            typename std::map<ProcessRank, BorderDistance>::const_iterator endIt = foreignIndexOverlap.end();
            for (int j = 0; it != endIt; ++it, ++j)
            {
                (*peersSendBuff_[peerRank][i])[2*j + 0] = it->first;
                (*peersSendBuff_[peerRank][i])[2*j + 1] = it->second;
            }
        }

        // send all messages
        indicesSendBuff_[peerRank]->send(peerRank);
        overlapIt = foreignOverlap.begin();
        for (int i = 0; overlapIt != overlapEndIt; ++overlapIt, ++i) {
            peersSendBuff_[peerRank][i]->send(peerRank);
        }
#endif // HAVE_MPI
    }

    void waitSendIndices_(int peerRank)
    {
        numIndicesSendBuff_[peerRank]->wait();
        numIndicesSendBuff_[peerRank].template reset<MpiBuffer<int> >(0);

        indicesSendBuff_[peerRank]->wait();
        indicesSendBuff_[peerRank].template reset<MpiBuffer<IndexDistanceNpeers> >(0);

        const ForeignOverlapWithPeer &foreignPeerOverlap
            = foreignOverlap_.foreignOverlapWithPeer(peerRank);
        ForeignOverlapWithPeer::const_iterator overlapIt = foreignPeerOverlap.begin();
        ForeignOverlapWithPeer::const_iterator overlapEndIt = foreignPeerOverlap.end();
        for (int i = 0; overlapIt != overlapEndIt; ++overlapIt, ++i) {
            peersSendBuff_[peerRank][i]->wait();
            peersSendBuff_[peerRank][i].template reset<MpiBuffer<int> >(0);
        }
    }

    void receiveIndicesFromPeer_(int peerRank)
    {
#if HAVE_MPI
        // receive the number of additional indices
        int numIndices = -1;
        MpiBuffer<int> numIndicesRecvBuff(1);
        numIndicesRecvBuff.receive(peerRank);
        numIndices = numIndicesRecvBuff[0];

        // receive the additional indices themselfs
        MpiBuffer<IndexDistanceNpeers> recvBuff(numIndices);
        recvBuff.receive(peerRank);
        for (int i = 0; i < numIndices; ++i) {
            int globalIdx = std::get<0>(recvBuff[i]);
            int borderDistance = std::get<1>(recvBuff[i]);
            int numPeers = std::get<2>(recvBuff[i]);

            int domesticIdx;
            if (borderDistance > 0) {
                // if the index is not on the border, add it to the
                // domestic overlap
                if (!globalIndices_.hasGlobalIndex(globalIdx)) {
                    // create and add a new domestic index
                    domesticIdx = globalIndices_.numDomestic();
                    borderDistance_.resize(domesticIdx + 1, std::numeric_limits<int>::max());
                    domesticOverlapByIndex_.resize(domesticIdx + 1);

                    globalIndices_.addIndex(domesticIdx, globalIdx);
                    domesticOverlapByIndex_[domesticIdx][peerRank] = borderDistance;
                    domesticOverlapWithPeer_[peerRank].insert(domesticIdx);
                }
                else {
                    domesticIdx = globalIndices_.globalToDomestic(globalIdx);
                }
            }
            else
                // border index
                domesticIdx = globalIndices_.globalToDomestic(globalIdx);

            borderDistance_[domesticIdx] = std::min(borderDistance, borderDistance_[domesticIdx]);

            // receive the peer ranks which see this index
            MpiBuffer<ProcessRank> peerRanksRecvBuff(2*numPeers);
            peerRanksRecvBuff.receive(peerRank);
            for (int j = 0; j < numPeers; ++j) {
                int seenBy = peerRanksRecvBuff[2*j + 0];
                int borderDistance2 = peerRanksRecvBuff[2*j + 1];
                if (seenBy != myRank_) {
                    domesticOverlapByIndex_[domesticIdx][seenBy] = borderDistance2;
                    domesticOverlapWithPeer_[seenBy].insert(domesticIdx);
                    peerSet_.insert(seenBy);
                }
            }
        }
#endif // HAVE_MPI
    }

    int myRank_;
    ForeignOverlap foreignOverlap_;

    DomesticOverlapByRank domesticOverlapWithPeer_;
    DomesticOverlapByIndex domesticOverlapByIndex_;
    std::vector<int> borderDistance_;

    std::map<ProcessRank, std::shared_ptr<MpiBuffer<int> > > numIndicesSendBuff_;
    std::map<ProcessRank, std::shared_ptr<MpiBuffer<IndexDistanceNpeers> > > indicesSendBuff_;
    std::map<ProcessRank, std::vector<std::shared_ptr<MpiBuffer<int> > > > peersSendBuff_;
    GlobalIndices globalIndices_;
    PeerSet peerSet_;
};

} // namespace Dumux

#endif
