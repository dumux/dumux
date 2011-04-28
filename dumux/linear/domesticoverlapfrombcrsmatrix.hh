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

    typedef std::pair<LocalIndex, BorderDistance> IndexDistance;
    typedef std::list<IndexDistance> DomesticOverlapWithPeer;
    typedef std::map<ProcessRank, DomesticOverlapWithPeer> DomesticOverlapByRank;
    typedef std::vector<std::map<ProcessRank, BorderDistance> > DomesticOverlapByIndex;

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
     * \brief Returns the number of processes which "see" a given
     *        index.
     */
    int numPeers(int domesticIdx) const
    { return numPeers_[domesticIdx]; };

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
    { return foreignOverlap_.peerSet(); }
    
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
        return foreignOverlap_.iAmMasterOf(domesticIdx);
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
        // resize the array which stores the number of peers for
        // each entry.
        numPeers_.resize(numDomestic(), 0);
        // for all local indices copy the number of processes from the
        // foreign overlap
        for (int i = 0; i < numLocal(); ++i) {
            numPeers_[i] = foreignOverlap_.numPeers(i);
        }

        // receive our overlap from the processes with higher rank
        PeerSet::const_iterator peerIt = peerSet().begin();
        PeerSet::const_iterator peerEndIt = peerSet().end();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            if (myRank_ < peerRank)
                receiveIndicesFromPeer_(peerRank);
        };

        // send the overlap indices to the peer processes with
        // lower rank
        peerIt = peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            if (myRank_ > peerRank)
                sendIndicesToPeer_(peerRank);
        };


        // receive our overlap from the processes with lower rank
        peerIt = peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            if (myRank_ > peerRank)
                receiveIndicesFromPeer_(peerRank);
        };

        // send the overlap indices to the peer processes with
        // higher rank
        peerIt = peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            if (myRank_ < peerRank)
                sendIndicesToPeer_(peerRank);
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

            int sendBuff[3] = {
                globalIndices_.domesticToGlobal(localIdx),
                borderDistance,
                0
            };
            if (foreignOverlap_.isBorderFor(peerRank, localIdx) &&
                peerRank < myRank_)
            {
                sendBuff[2] = foreignOverlap_.numNonBorderPeers(localIdx);
            }
            else
                sendBuff[2] = numPeers_[localIdx];
            
            MPI_Bsend(sendBuff, // buff
                      3, // count
                      MPI_INT, // data type
                      peerRank, 
                      0, // tag
                      MPI_COMM_WORLD); // communicator
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
                    numPeers_.push_back(0);
                }
                else {
                    domesticIdx = globalIndices_.globalToDomestic(globalIdx);
                }
                
                domesticOverlapWithPeer_[peerRank].push_back(IndexDistance(domesticIdx, borderDistance));
            }
            else 
                // border index
                domesticIdx = globalIndices_.globalToDomestic(globalIdx);
            
            // calculate the number of peers
            if (isBorder(domesticIdx) && peerRank > myRank_)
                numPeers_[domesticIdx] += numPeers;
            else
                numPeers_[domesticIdx] = numPeers;
        }
    }

    int myRank_;
    ForeignOverlap foreignOverlap_;

    DomesticOverlapByRank domesticOverlapWithPeer_;
    std::vector<int> numPeers_;

    GlobalIndices globalIndices_;
};

} // namespace Dumux

#endif
