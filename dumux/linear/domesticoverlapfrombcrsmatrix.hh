// $Id$
/*****************************************************************************
 *   Copyright (C) 2011 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
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

#include "foreignoverlapfrombcrsmatrix.hh"
#include "globalindices.hh"

#include <algorithm>
#include <list>
#include <set>
#include <map>

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
    typedef Index PeerIndex;
    typedef Index LocalIndex;
    typedef std::pair<LocalIndex, ProcessRank> IndexRank;
    typedef std::tuple<Index, BorderDistance, int> IndexDistanceNpeers;
    typedef std::list<IndexRank> SeedList;

    typedef std::set<ProcessRank> PeerSet;

    typedef std::list<IndexDistanceNpeers> ForeignOverlapWithPeer;
    typedef std::vector<std::map<ProcessRank, BorderDistance> > ForeignOverlapByIndex;

    typedef std::pair<LocalIndex, BorderDistance> IndexDistance;
    typedef std::set<LocalIndex> DomesticOverlapWithPeer;
    typedef std::map<ProcessRank, DomesticOverlapWithPeer> DomesticOverlapByRank;
    typedef std::vector<std::set<ProcessRank> > DomesticOverlapByIndex;

    typedef std::tuple<LocalIndex, PeerIndex, ProcessRank> LocalindexPeerindexPeerrank;
    typedef std::list<LocalindexPeerindexPeerrank> BorderList;

    /*!
     * \brief Constructs the foreign overlap given a BCRS matrix and
     *        an initial list of border indices.
     */
    DomesticOverlapFromBCRSMatrix(const BCRSMatrix &A,
                                  const BorderList &borderList,
                                  int overlapSize)
        : foreignOverlap_(A, borderList, overlapSize)
        , globalIndices_(foreignOverlap_)
    {
        myRank_ = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &myRank_);

        buildDomesticOverlap_(A);
    }

    /*!
     * \brief Return the border list.
     *
     * The border list is the list of (local index, peer index, peer
     * rank) triples for all indices on a process border.
     */
    const BorderList& borderList() const
    { return foreignOverlap_.borderList(); };
     

    /*!
     * \brief Returns true iff a domestic index is a border index.
     */
    bool isBorder(int domesticIdx) const
    { 
        return isLocal(domesticIdx) && foreignOverlap_.isBorder(domesticIdx);
    };

    /*!
     * \brief Returns true iff a domestic index is a front index.
     */
    bool isFront(int domesticIdx) const
    { return borderDistance_[domesticIdx] == foreignOverlap_.overlapSize(); };

    /*!
     * \brief Returns the number of processes which "see" a given
     *        index.
     */
    int numPeers(int domesticIdx) const
    { return domesticOverlapByIndex_[domesticIdx].size(); };

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
    { return foreignOverlap_.numLocal(); };

    /*!
     * \brief Returns the number domestic indices.
     *
     * The domestic indices are defined as the process' local indices
     * plus its copies of indices in the overlap regions
     */
    int numDomestic() const
    { return globalIndices_.numDomestic(); };

    /*!
     * \brief Return true if a domestic index is local for the process
     *        (i.e. interior or border)
     */
    bool isLocal(int domesticIdx) const
    { return domesticIdx < numLocal(); };

    /*!
     * \brief Return true iff the current process is the master of a
     *        given domestic index.
     */
    bool iAmMasterOf(int domesticIdx) const
    { 
        if (!isLocal(domesticIdx))
            return false;
        return foreignOverlap_.masterOf(domesticIdx) == myRank_;
    };

    /*!
     * \brief Print the foreign overlap for debugging purposes.
     */
    void print() const
    {
        globalIndices_.print();
    };

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
        for (int i = 0; i < numLocal(); ++i) {
            const std::map<ProcessRank, BorderDistance> &idxOverlap =
                foreignOverlap_.foreignOverlapByIndex(i);
            std::map<ProcessRank, BorderDistance>::const_iterator it = idxOverlap.begin();
            std::map<ProcessRank, BorderDistance>::const_iterator endIt = idxOverlap.end();
            for (; it != endIt; ++it) {
                domesticOverlapByIndex_[i].insert(it->first);
                domesticOverlapWithPeer_[it->first].insert(i);
            }
        }

        PeerSet::const_iterator peerIt;
        PeerSet::const_iterator peerEndIt = foreignOverlap_.peerSet().end();

        // send the overlap indices to all peer processes 
        peerIt = foreignOverlap_.peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            sendIndicesToPeer_(peerRank);
        };

        // receive our overlap from the processes to all peer processes
        peerIt = foreignOverlap_.peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            receiveIndicesFromPeer_(peerRank);
        };
    };

    void sendIndicesToPeer_(int peerRank)
    {          
        const ForeignOverlapWithPeer &foreignOverlap
            = foreignOverlap_.foreignOverlapWithPeer(peerRank);

        // first, send a message containing the number of additional
        // indices stemming from the overlap (i.e. without the border
        // indices)
        int numIndices = foreignOverlap.size();
        MPI_Send(&numIndices, // buff
                 1, // count
                 MPI_INT, // data type
                 peerRank, 
                 0, // tag
                 MPI_COMM_WORLD); // communicator

        // then send the additional indices themselfs
        ForeignOverlapWithPeer::const_iterator overlapIt = foreignOverlap.begin();
        ForeignOverlapWithPeer::const_iterator overlapEndIt = foreignOverlap.end();
        for (; overlapIt != overlapEndIt; ++overlapIt) {
            int localIdx = std::get<0>(*overlapIt);
            int borderDistance = std::get<1>(*overlapIt);

            const std::map<ProcessRank, BorderDistance> &foreignIndexOverlap
                = foreignOverlap_.foreignOverlapByIndex(localIdx);

            int numPeers = foreignIndexOverlap.size();

            int sendBuff[3] = {
                globalIndices_.domesticToGlobal(localIdx),
                borderDistance,
                numPeers
            };

            MPI_Bsend(sendBuff, // buff
                      3, // count
                      MPI_INT, // data type
                      peerRank, 
                      0, // tag
                      MPI_COMM_WORLD); // communicator

            int *peerRanks = new int[numPeers];
            typename std::map<ProcessRank, BorderDistance>::const_iterator it = foreignIndexOverlap.begin();
            typename std::map<ProcessRank, BorderDistance>::const_iterator endIt = foreignIndexOverlap.end();
            for (int i = 0; it != endIt; ++it, ++i)
            { peerRanks[i] = it->first; };
            
            MPI_Bsend(peerRanks, // buff
                      numPeers, // count
                      MPI_INT, // data type
                      peerRank, 
                      0, // tag
                      MPI_COMM_WORLD); // communicator
            delete[] peerRanks;
        };
    }

    void receiveIndicesFromPeer_(int peerRank)
    {
        // receive the number of additional indices
        int numIndices = -1;
        MPI_Recv(&numIndices, // buff
                 1, // count
                 MPI_INT, // data type
                 peerRank, 
                 0, // tag
                 MPI_COMM_WORLD, // communicator
                 MPI_STATUS_IGNORE);
        
        // receive the additional indices themselfs
        for (int i = 0; i < numIndices; ++i) {
            int recvBuff[3];
            MPI_Recv(recvBuff, // buff
                     3, // count
                     MPI_INT, // data type
                     peerRank, 
                     0, // tag
                     MPI_COMM_WORLD, // communicator
                     MPI_STATUS_IGNORE);

            int globalIdx = recvBuff[0];
            int borderDistance = recvBuff[1];
            int numPeers = recvBuff[2];
            int domesticIdx;
            if (borderDistance > 0) {
                // if the index is not on the border, add it to the
                // domestic overlap
                if (!globalIndices_.hasGlobalIndex(globalIdx)) {
                    // create and add a new domestic index
                    domesticIdx = globalIndices_.numDomestic();
                    globalIndices_.addIndex(domesticIdx, globalIdx);
                    
                    int nDom = globalIndices_.numDomestic();
                    domesticOverlapByIndex_.resize(nDom);
                    borderDistance_.resize(nDom, std::numeric_limits<int>::max());
                }
                else {
                    domesticIdx = globalIndices_.globalToDomestic(globalIdx);
                }
            }
            else 
                // border index
                domesticIdx = globalIndices_.globalToDomestic(globalIdx);

            borderDistance_[domesticIdx] = std::min(borderDistance, borderDistance_[domesticIdx]);

            int *peerRanks = new int[numPeers];
            MPI_Recv(peerRanks, // buff
                     numPeers, // count
                     MPI_INT, // data type
                     peerRank, 
                     0, // tag
                     MPI_COMM_WORLD, // communicator
                     MPI_STATUS_IGNORE);

            domesticOverlapByIndex_[domesticIdx].insert(peerRank);
            domesticOverlapWithPeer_[peerRank].insert(domesticIdx);
            for (int i = 0; i < numPeers; ++i) {
                if (peerRanks[i] != myRank_) {
                    domesticOverlapByIndex_[domesticIdx].insert(peerRanks[i]);
                    domesticOverlapWithPeer_[peerRanks[i]].insert(domesticIdx);
                    peerSet_.insert(peerRanks[i]);
                }
            }
            delete [] peerRanks;
        }
    }

    int myRank_;
    ForeignOverlap foreignOverlap_;

    DomesticOverlapByRank domesticOverlapWithPeer_;
    DomesticOverlapByIndex domesticOverlapByIndex_;
    std::vector<int> borderDistance_;

    GlobalIndices globalIndices_;
    PeerSet peerSet_;
};

} // namespace Dumux

#endif
