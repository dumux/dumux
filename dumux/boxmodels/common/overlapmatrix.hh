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
 * \brief A BCRS matrix which creates an algebraic overlap of
 *        arbitrary size.
 */
#ifndef DUMUX_OVERLAP_MATRIX_HH
#define DUMUX_OVERLAP_MATRIX_HH

#include <dune/grid/common/datahandleif.hh>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <algorithm>
#include <list>
#include <set>
#include <map>


namespace Dumux {

/*!
 * \brief Uses communication on the grid to find the initial seed list
 *        of indices.
 *
 * \todo implement this class generically. For this, it must be
 *       possible to query the mapper whether it contains entities of
 *       a given codimension without the need to hand it an actual
 *       entity.
 */
template <class GridView, class VertexMapper>
class VertexBorderListFromGrid : public Dune::CommDataHandleIF<VertexBorderListFromGrid<GridView, VertexMapper>,
                                                               int >
{
    typedef int ProcessRank;
    typedef int Index;
    typedef Index LocalIndex;
    typedef Index PeerIndex;
    typedef std::tuple<LocalIndex, PeerIndex, ProcessRank> LindexPindexRank;
    typedef std::list<LindexPindexRank> BorderList;

public:
    VertexBorderListFromGrid(const GridView &gv, 
                             const VertexMapper &map)
        : gv_(gv), map_(map)
    {
        gv.communicate(*this, 
                       Dune::InteriorBorder_InteriorBorder_Interface, 
                       Dune::ForwardCommunication);
    };

    // data handle methods
    bool contains (int dim, int codim) const
    { return dim == codim; }
    
    bool fixedsize(int dim, int codim) const
    { return true; }
    
    template<class EntityType> 
    size_t size(const EntityType &e) const
    { return 2; };

    template<class MessageBufferImp, class EntityType> 
    void gather(MessageBufferImp &buff, const EntityType &e) const 
    { 
        buff.write(gv_.comm().rank()); 
        buff.write(map_.map(e)); 
    };

    template<class MessageBufferImp, class EntityType> 
    void scatter(MessageBufferImp &buff, const EntityType &e, size_t n)
    {
        int peerRank;
        int peerVertIdx;

        buff.read(peerRank);
        buff.read(peerVertIdx);

        int localVertIdx = map_.map(e);
        borderList_.push_back(LindexPindexRank(localVertIdx, peerVertIdx, peerRank));
    };
    
    // Access to the initial seed list.
    const BorderList &borderList() const
    { return borderList_; }

private:
    const GridView gv_;
    const VertexMapper &map_;
    BorderList borderList_;
};


/*!
 * \brief Creates a overlap given an arbitrary BCRS matrix and an
 *        initial list of rows on the process border.
 */
template<class BCRSMatrix>
class OverlapFromBCRSMatrix
{
    OverlapFromBCRSMatrix(const OverlapFromBCRSMatrix &A)
    {}

public:
    typedef int ProcessRank;
    typedef int BorderDistance;
    typedef int Index;
    typedef Index PeerIndex;
    typedef Index DomesticIndex;
    typedef std::pair<DomesticIndex, ProcessRank> IndexRank;
    typedef std::pair<PeerIndex, BorderDistance> IndexDistance;
    typedef std::tuple<Index, BorderDistance, int>  IndexDistanceNpeers;
    typedef std::list<IndexRank> SeedList;
    typedef std::pair<ProcessRank, int> ProcessDistance;

    typedef std::set<ProcessRank> PeerSet;
    typedef std::list<IndexDistanceNpeers> ForeignOverlapWithPeer;
    typedef std::map<ProcessRank,  ForeignOverlapWithPeer> ForeignOverlapByRank;
    typedef std::vector<std::map<ProcessRank, BorderDistance> > ForeignOverlapByIndex;

    typedef std::vector<IndexDistanceNpeers> DomesticOverlapWithPeer;
    typedef std::map<ProcessRank,  DomesticOverlapWithPeer> DomesticOverlap;
    typedef std::map<ProcessRank, std::set<PeerIndex> > DomesticFront;
    typedef std::map<ProcessRank, std::set<DomesticIndex> > ForeignFront;
    typedef std::map<PeerIndex, DomesticIndex> BorderWithPeer;
    typedef std::map<ProcessRank, BorderWithPeer> DomesticBorder;

    typedef std::tuple<DomesticIndex, PeerIndex, ProcessRank> LindexPindexRank;
    typedef std::list<LindexPindexRank> BorderList;

    OverlapFromBCRSMatrix(const BCRSMatrix &A,
                          const BorderList &borderList,
                          int overlapSize)
    {
        // calculate the border list. From this, create an initial
        // seed list of indices which are in the overlap.
        SeedList initialSeedList;
        buildBorder_(initialSeedList, borderList);

        // find the set of processes which have an overlap with the
        // local processes. (i.e. the set of processes which we will
        // have to communicate to.)
        buildPeerSet_(initialSeedList);

        // calculate the foreign overlap for the local partition,
        // i.e. find the distance of each row from the seed set.
        foreignOverlapByIndex_.resize(A.N());
        extendForeignOverlap_(A, initialSeedList, 0, overlapSize);

        // group foreign overlap by peer process rank
        groupForeignOverlapByRank_(); 
        
        // calculate the domestic overlap (i.e. all overlap indices in
        // foreign processes which the current process overlaps.)
        // This requires communication via MPI.
        buildDomesticOverlap_();

        // calculate the border and front lists
        buildFront_(overlapSize);

        // set the number of local indices to the number of rows of
        // the BCRS matrix
        numLocalIndices_ = A.N();

        //  fill the (peerRank, peerIndex) -> domestic index map
        buildForeignToDomesticMap_();

        // fills the array which stores the number of peers for each
        // domestic index.
        buildNumOverlapPeers_();
    }

    int numLocal() const
    { return numLocalIndices_; }

    int numDomestic() const
    { return numDomesticIndices_; }

    bool isLocal(int idx) const
    { return idx < numLocalIndices_; };

    bool isBorder(int peerRank, int foreignIdx) const
    { return domesticBorder_.find(peerRank)->second.count(foreignIdx) > 0; };

    bool isFront(int peerRank, int foreignIdx) const
    { return domesticFront_.find(peerRank)->second.count(foreignIdx) > 0; };

    bool isForeignFront(int peerRank, int domesticIdx) const
    { return foreignFront_.find(peerRank)->second.count(domesticIdx) > 0; };

    const ForeignOverlapWithPeer &foreignOverlapWithPeer(int peerRank) const
    { return foreignOverlapByRank_.find(peerRank)->second; }

    const DomesticOverlapWithPeer &domesticOverlapWithPeer(int peerRank) const
    { return domesticOverlap_.find(peerRank)->second; }

    const BorderWithPeer &borderWithPeer(int peerRank) const
    { return domesticBorder_.find(peerRank)->second; }

    const PeerSet &peerSet() const
    { return peerSet_; }
    
    int foreignToDomesticIndex(int peerRank, int foreignIdx) const
    { 
        assert(foreignToDomesticMap_.count(peerRank) > 0);
        //assert(foreignToDomesticMap_[peerRank].count(foreignIdx) > 0);

        return foreignToDomesticMap_.find(peerRank)->second.find(foreignIdx)->second;
    };

    int domesticToForeignIndex(int domesticIdx, int peerRank) const
    { 
        assert(domesticOffset_.count(peerRank) > 0);
        //assert(domesticIdx >= numLocalIndices_ + domesticOffset_[peerRank]);
        //assert(domesticIdx < numLocalIndices_ + domesticOffset_[peerRank] + domesticOverlap_[peerRank].size());

        return domesticToForeignMap_.find(peerRank)->second.find(domesticIdx)->second;
    };

    int numOverlapPeers(int domesticIdx) const
    { return numOverlapPeers_[domesticIdx]; };

    void printOverlap() const
    {
        int myRank;
        MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
        std::cout << "Number of ranks with overlap: " << peerSet_.size() << "\n";

        DomesticOverlap::const_iterator it = domesticOverlap_.begin();
        DomesticOverlap::const_iterator endIt = domesticOverlap_.end();
        for (; it != endIt; ++it) {
            int peerRank = it->first;
            std::cout << "Rank " << myRank << " overlaps with rank " << peerRank << " at: ";

            std::vector<IndexDistanceNpeers>::const_iterator rowIt = it->second.begin();
            std::vector<IndexDistanceNpeers>::const_iterator rowEndIt = it->second.end();
            for (; rowIt != rowEndIt; ++rowIt) {
                int peerIdx = std::get<0>(*rowIt);
                int domIdx = foreignToDomesticIndex(peerRank, peerIdx);
                int borderDist = std::get<1>(*rowIt);
                int numRanks = std::get<2>(*rowIt);

                int peerIdx2 = domesticToForeignIndex(domIdx, peerRank);
                assert(isBorder(peerRank, peerIdx) || peerIdx2 == peerIdx);

                std::cout <<  peerIdx << "(dist=" << borderDist 
                          << ",np=" << numRanks
                          << ",domIdx=" << domIdx
                          << ",pi2="<<peerIdx2
                          << ") ";
            };
            std::cout << "\n";
        }

/*
        ForeignOverlapByRank::const_iterator it = foreignOverlapByRank_.begin();
        ForeignOverlapByRank::const_iterator endIt = foreignOverlapByRank_.end();
        for (; it != endIt; ++it) {
            std::cout << "Overlap rows(distance) for rank " << it->first << ": ";
            
            std::list<IndexDistance>::const_iterator rowIt = it->second.begin();
            std::list<IndexDistance>::const_iterator rowEndIt = it->second.end();
            for (; rowIt != rowEndIt; ++rowIt) {
                std::cout << rowIt->first << "(" << rowIt->second << ") ";
            };
            std::cout << "\n";
        }
*/
/*
        int n = foreignOverlapByIndex_.size();
        for (int i = 0; i < n; ++i) {
            const std::map<int, int> &rowSet = foreignOverlapByIndex_[i];
            if (rowSet.size() == 0)
                continue;
            
            std::cout << "Overlap ranks(distance) for row " << i << ": ";
            
            std::map<int, int>::const_iterator it = rowSet.begin();
            std::map<int, int>::const_iterator endIt = rowSet.end();
            for (; it != endIt; ++it) {
                std::cout << it->first << "(" << it->second << ")" << " ";
            };
            std::cout << "\n";
        };
*/
    };
                  
protected:
    // Calculate the set of peer processes from the initial seed list.
    void buildPeerSet_(const SeedList &seedList)
    {
        SeedList::const_iterator it = seedList.begin();
        SeedList::const_iterator endIt = seedList.end();
        for (; it != endIt; ++it)
            peerSet_.insert(it->second);
    }
    
    // extend the foreign overlaps by one level. this uses a greedy
    // algorithm.
    void extendForeignOverlap_(const BCRSMatrix &A,
                               SeedList &seedList,
                               int overlapDistance,
                               int overlapSize)
    {
        // add all processes in the seed rows of the current overlap
        // level
        SeedList::const_iterator it = seedList.begin();
        SeedList::const_iterator endIt = seedList.end();
        for (; it != endIt; ++it) {
            if (foreignOverlapByIndex_[it->first].count(it->second) == 0) {
                foreignOverlapByIndex_[it->first][it->second] = overlapDistance;
            }
        }

        // if we have an overlap of 0, we're finished and break the
        // recursion
        if (overlapSize == overlapDistance)
            return;
            
        // find the seed list for the next overlap level using the
        // seed set for the current level
        SeedList nextSeedList;
        it = seedList.begin();
        for (; it != endIt; ++it) {
            int rowIdx = it->first;
            int processIdx = it->second;
                
            // find all column indices in the row. The indices of the
            // columns are the additional indices of the overlap which
            // we would like to add
            typedef typename BCRSMatrix::ConstColIterator ColIterator;
            ColIterator colIt = A[rowIdx].begin();
            ColIterator colEndIt = A[rowIdx].end();
            for (; colIt != colEndIt; ++colIt) {
                int newIdx = colIt.index();
                    
                // if the process is already is in the overlap of the
                // column index, ignore this column index!
                if (foreignOverlapByIndex_[newIdx].count(processIdx) > 0)
                    continue;
                    
                // check whether the new index is already in the next seed list
                IndexRank newPair(newIdx, processIdx);
                if (std::find(nextSeedList.begin(), 
                              nextSeedList.end(),
                              newPair) != nextSeedList.end())
                    continue; // we already have this pair                    

                // add the current processes to the seed list for the
                // next overlap level
                nextSeedList.push_back(newPair);
            }
        }
        
        // clear the old seed list to save some memory
        seedList.clear();
        
        // Perform the same excercise for the next overlap level
        extendForeignOverlap_(A,
                              nextSeedList, 
                              overlapDistance + 1, 
                              overlapSize);
    };

    // assuming that the foreign overlap has been created for each
    // local index, this method groups the foreign overlap by peer
    // process rank
    void groupForeignOverlapByRank_()
    {
        int nIndices = foreignOverlapByIndex_.size();
        for (int i = 0; i < nIndices; ++i)
        {
            std::map<ProcessRank, BorderDistance>::const_iterator it = 
                foreignOverlapByIndex_[i].begin();
            std::map<ProcessRank, BorderDistance>::const_iterator endIt = 
                foreignOverlapByIndex_[i].end();
            int nRanks = foreignOverlapByIndex_[i].size();
            for (; it != endIt; ++it)  {
                IndexDistanceNpeers tmp(i, it->second, nRanks);
                foreignOverlapByRank_[it->first].push_back(tmp);
            };
        };
    }
    
    // build the domestic overlap of the process, i.e. the overlapping
    // indices of the current process within the peer processes
    void buildDomesticOverlap_()
    {
        ///////
        // Asyncronously send all indices to the peer ranks
        ///////
        int myRank;
        MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

        MPI_Request sendRequests[peerSet_.size()][2];

        int *sendMsgBuff[peerSet_.size()];
        std::set<int>::iterator it = peerSet_.begin();
        std::set<int>::iterator endIt = peerSet_.end();
        for (int i = 0; it != endIt; ++it, ++i) {
            int peerRank = *it;
            assert(peerRank != myRank);

            int sendMsgSize;
            sendMsgBuff[i] = packSendMessage_(peerRank, sendMsgSize);
            
            // send the size of the overlap message buffer to the
            // peer rank
            MPI_Isend(&sendMsgSize, // pointer to user data 
                      1, // size of user data array
                      MPI_INT, // type of user data
                      peerRank,  // peer rank
                      0, // identifier
                      MPI_COMM_WORLD, // communicator
                      &sendRequests[i][0]); // request object
            // send the overlap indices and distances to the peer
            // rank
            MPI_Isend(sendMsgBuff[i], // pointer to user data 
                      sendMsgSize, // size of user data array
                      MPI_INT, // type of user data
                      peerRank,  // peer rank
                      0, // identifier
                      MPI_COMM_WORLD, // communicator
                      &sendRequests[i][1]); // request object
        };

        ///////
        // Receive all peer indices using syncronous receives
        ///////
        it = peerSet_.begin();
        for (; it != endIt; ++it) {
            assert(*it != myRank);
            
            // receive size of overlap message buffer
            int recvMsgSize = -1;
            MPI_Recv(&recvMsgSize, 1, MPI_INT, *it, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            // retrieve the actual overlap array message
            int *recvMsgBuff = new int[recvMsgSize];
            MPI_Recv(recvMsgBuff, recvMsgSize, MPI_INT, *it, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            // unpack the overlap message
            unpackRecvMessage_(*it, recvMsgBuff, recvMsgSize);
            delete[] recvMsgBuff;
        };

        ///////
        // Wait for all send requests to complete, (this should not be
        // necessary, since we used syncronous receives, but it might
        // avoid memory leaks inside the MPI library) and delete the
        // send buffers
        ///////
        it = peerSet_.begin();
        for (int i = 0; it != endIt; ++it, ++i) {
            MPI_Wait(&sendRequests[i][0], MPI_STATUS_IGNORE);
            MPI_Wait(&sendRequests[i][1], MPI_STATUS_IGNORE);

            // free the memory used by the message buffer
            delete[] sendMsgBuff[i];
        }
    }

    // pack the overlap message buffer send to a peer
    int *packSendMessage_(int peerRank, int &msgSize) 
    {
        int n = foreignOverlapByRank_[peerRank].size();
        msgSize = 3*n;
        int *buff = new int[msgSize];

        std::list<IndexDistanceNpeers>::const_iterator it = foreignOverlapByRank_[peerRank].begin();
        std::list<IndexDistanceNpeers>::const_iterator endIt = foreignOverlapByRank_[peerRank].end();
        for (int i = 0; i < n; ++i, ++it) {
            assert(std::get<0>(*it) >= 0);
            assert(std::get<1>(*it) >= 0);
            buff[3*i + 0] = std::get<0>(*it); // row index
            buff[3*i + 1] = std::get<1>(*it); // border distance
            buff[3*i + 2] = std::get<2>(*it); // number of ranks

            assert(buff[3*i + 0] >= 0);
            assert(buff[3*i + 1] >= 0);
            assert(buff[3*i + 2] > 0);
        };

        return buff;
    };

    // unpack the overlap message buffer received from a peer
    void unpackRecvMessage_(int peerRank, int *overlapMsg, int &msgSize) 
    {
        int n = msgSize / 3;
        
        // create the domestic overlap for the peer rank
        for (int i = 0; i < n; ++i) {
            int peerIdx = overlapMsg[3*i + 0];
            int borderDistance = overlapMsg[3*i + 1];
            int numRanks = overlapMsg[3*i + 2];
            assert(peerIdx >= 0);
            assert(borderDistance >= 0);
            assert(numRanks > 0);
            IndexDistanceNpeers tmp(peerIdx, borderDistance, numRanks);
            domesticOverlap_[peerRank].push_back(tmp);
        };
    };

    // calculate the domestic border indices given the initial seed list
    void buildBorder_(SeedList &initialSeedList, const BorderList &borderList)
    {
        BorderList::const_iterator it = borderList.begin();
        BorderList::const_iterator endIt = borderList.end();
        for (; it != endIt; ++it) {
            int localIdx = std::get<0>(*it);
            int peerIdx = std::get<1>(*it);
            int peerRank = std::get<2>(*it);
            
            initialSeedList.push_back(IndexRank(localIdx, peerRank));
            domesticBorder_[peerRank][peerIdx] = localIdx;
        };
    };
    
    // assuming that the domestic overlap has already been calulated,
    // calculate the front indices from the overlap distance
    void buildFront_(int overlapSize)
    {
        int myRank;
        MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

        std::set<int>::iterator it = peerSet_.begin();
        std::set<int>::iterator endIt = peerSet_.end();
        for (; it != endIt; ++it) {
            int peerRank = *it;
            assert(peerRank != myRank);
            
            // find the our own front (i.e. the dometic front)
            std::vector<IndexDistanceNpeers >::iterator oit = domesticOverlap_[peerRank].begin();
            std::vector<IndexDistanceNpeers >::iterator oEndIt = domesticOverlap_[peerRank].end();
            for (; oit != oEndIt; ++oit) {
                int peerIdx = std::get<0>(*oit);
                int peerDistance = std::get<1>(*oit);
                
                if (peerDistance == overlapSize) 
                    // front index. TODO: this is only true for
                    // indices which are not on the boundary of the
                    // peer proces' domain. we probably have to
                    // communicate to get the real front set!
                    domesticFront_[peerRank].insert(peerIdx);
            }

            // find the front of the peers (i.e. the foreign front)
            std::list<IndexDistanceNpeers >::iterator foit = foreignOverlapByRank_[peerRank].begin();
            std::list<IndexDistanceNpeers >::iterator foEndIt = foreignOverlapByRank_[peerRank].end();
            for (; foit != foEndIt; ++foit) {
                int domesticIdx = std::get<0>(*foit);
                int peerDistance = std::get<1>(*foit);
                
                if (peerDistance == overlapSize) 
                    // front index. TODO: this is only true for
                    // indices which are not on the boundary of the
                    // peer proces' domain. we probably have to
                    // communicate to get the real front set!
                    foreignFront_[peerRank].insert(domesticIdx);
            }
        };
    };
    
    // builds the map to translate from a (peer rank, foreign index)
    // pair to the domestic index
    void buildForeignToDomesticMap_()
    {
        // calculate the offsets for the additional domestic indices
        numDomesticIndices_ = 0;
        DomesticOverlap::iterator ovrlpIt = domesticOverlap_.begin();
        DomesticOverlap::iterator ovrlpEndIt = domesticOverlap_.end();
        for (; ovrlpIt != ovrlpEndIt; ++ovrlpIt) {
            int peerRank = ovrlpIt->first;

            // calculate the number of border indices for the current
            // peer rank.
            int numBorder = 0;
            int numIndices = ovrlpIt->second.size();
            for (int i = 0; i < numIndices - numBorder; ++i) {
                int borderDist = std::get<1>(ovrlpIt->second[i]);

                if (borderDist == 0) {
                    // move this index to the back
                    std::swap(ovrlpIt->second[i],
                              ovrlpIt->second[numIndices - 1 - numBorder]);
                    
                    // increment border count
                    ++ numBorder;
                    --i; // but don't look at the next entry in the
                         // next iteration
                    continue;
                }
            }
                        
            //borderIndices_[peerRank] = nBorder;
            domesticOffset_[peerRank] = numDomesticIndices_;
            numDomesticIndices_ += numIndices - numBorder;
        }
        numDomesticIndices_ += numLocalIndices_;

        // add the border rows to the map
        DomesticBorder::const_iterator brdrIt = domesticBorder_.begin();
        DomesticBorder::const_iterator brdrEndIt = domesticBorder_.end();
        for (; brdrIt != brdrEndIt; ++brdrIt) {
            int peerRank = brdrIt->first;
            BorderWithPeer::const_iterator it = brdrIt->second.begin();
            BorderWithPeer::const_iterator endIt = brdrIt->second.end();
            for (; it != endIt; ++it) {
                int peerIdx = it->first;
                int localIdx = it->second;
                foreignToDomesticMap_[peerRank][peerIdx] = localIdx;
                domesticToForeignMap_[peerRank][localIdx] = peerIdx;
            }
        };        

        // add the additional indices to the {(foreignIndex, peerRank)
        // -> domestic index} map
        ovrlpIt = domesticOverlap_.begin();
        for (; ovrlpIt != ovrlpEndIt; ++ovrlpIt) {
            int peerRank = ovrlpIt->first;

            int numIndices = ovrlpIt->second.size();
            for (int i = 0, j = 0; i < numIndices; ++i) {
                int peerIdx = std::get<0>(ovrlpIt->second[i]);
                int borderDist = std::get<1>(ovrlpIt->second[i]);

                if (borderDist == 0) {
                    // ignore border indices since they are not
                    // additional degrees of freedom
                    continue;
                }

                int domIdx = numLocalIndices_ + domesticOffset_[peerRank] + j;
                foreignToDomesticMap_[peerRank][peerIdx] = domIdx;
                domesticToForeignMap_[peerRank][domIdx] = peerIdx;
                assert(domIdx >= numLocalIndices_);
                assert(domIdx < numDomesticIndices_);
                ++ j;
            };
        };
    };

    // calculates the number of peers which overlap any given domestic
    // index
    void buildNumOverlapPeers_()
    {
        int n = numDomestic();
        numOverlapPeers_.resize(n);
        std::fill(numOverlapPeers_.begin(),
                  numOverlapPeers_.end(),
                  0);

        // deal with the domestic overlap
        DomesticOverlap::iterator ovrlpIt = domesticOverlap_.begin();
        DomesticOverlap::iterator ovrlpEndIt = domesticOverlap_.end();
        for (; ovrlpIt != ovrlpEndIt; ++ovrlpIt) {
            int peerRank = ovrlpIt->first;
            
            const DomesticOverlapWithPeer &peerOverlap 
                = ovrlpIt->second;

            int numIndices = ovrlpIt->second.size();
            for (int i = 0; i < numIndices; ++i) {
                int rowIdx = std::get<0>(peerOverlap[i]);
                int numPeers = std::get<2>(peerOverlap[i]);
                
                int domRowIdx = foreignToDomesticIndex(peerRank, rowIdx);

                numOverlapPeers_[domRowIdx] = numPeers;
            };
        }

        // deal with the foreign overlap
        for (int i = 0; i < numLocalIndices_; ++i)
            numOverlapPeers_[i] = foreignOverlapByIndex_[i].size();
    };

    // set of processes with which we have to communicate
    PeerSet peerSet_;
    
    // stores the set of process ranks which are in the overlap for a
    // given row index "owned" by the current rank. The second value
    // store the distance from the nearest process border.
    ForeignOverlapByIndex foreignOverlapByIndex_;

    // stores a list of foreign overlap indices for each rank
    ForeignOverlapByRank foreignOverlapByRank_;
    
    // stores the overlap the current processes has in the domain
    // "owned" by some peer ranks
    DomesticOverlap domesticOverlap_;

    // stores the "front" indices of the domestic overlap
    DomesticFront domesticFront_;

    // stores the "front" indices of the foreign overlap
    ForeignFront foreignFront_;

    // stores the "border" indices of the domestic overlap
    DomesticBorder domesticBorder_;

    // number of local indices
    int numLocalIndices_;
    
    // number of domestic indices (i.e. local + additional indices for
    // the domestic overlap
    int numDomesticIndices_;
    
    // offset for the overlap indices for each peer process in number
    // of indices after the last "owned" index
    std::map<ProcessRank, int> domesticOffset_;

    // maps a peer local index to a domestic local index
    std::map<ProcessRank, std::map<PeerIndex, DomesticIndex> > foreignToDomesticMap_;
    std::map<ProcessRank, std::map<DomesticIndex, PeerIndex> > domesticToForeignMap_;

    // stores the number of peers which overlap any given domestic
    // index
    std::vector<int> numOverlapPeers_;
};

template <class BCRSMatrix, class Overlap>
class OverlapBCRSMatrix : public BCRSMatrix
{
    typedef BCRSMatrix ParentType;
    typedef typename Overlap::Index RowIndex;
    typedef typename Overlap::Index ColIndex;
    typedef typename Overlap::PeerSet PeerSet;
    typedef typename Overlap::ProcessRank ProcessRank;
    typedef typename Overlap::ForeignOverlapWithPeer ForeignOverlapWithPeer;
    typedef typename Overlap::DomesticOverlapWithPeer DomesticOverlapWithPeer;
    typedef std::vector<ColIndex> PeerColumns;
    typedef std::pair<RowIndex,  PeerColumns> PeerRow;
    typedef std::vector<PeerRow> PeerRows;
    typedef std::map<ProcessRank, PeerRows> OverlapStructure;

public:
    typedef typename ParentType::RowIterator RowIterator;
    typedef typename ParentType::ColIterator ColIterator;
    typedef typename ParentType::ConstColIterator ConstColIterator;
    typedef typename ParentType::field_type field_type;
    typedef typename ParentType::block_type block_type;

    OverlapBCRSMatrix(const BCRSMatrix &M, const Overlap &overlap)
        : overlap_(&overlap)
    {
        // build the overlapping matrix from the non-overlapping
        // matrix and the overlap
        build_(M);
    };
    
    OverlapBCRSMatrix(const OverlapBCRSMatrix &M)
        : ParentType(M), 
          overlap_(M.overlap_),
          numOverlapEntries_(M.numOverlapEntries_),
          overlapEntries_(M.overlapEntries_)
    {
    };

    void set(const BCRSMatrix &M)
    {
        // set everything to 0
        ParentType::operator=(0);

        // assign the local rows
        for (int rowIdx = 0; rowIdx < M.N(); ++rowIdx) {
            ConstColIterator colIt = M[rowIdx].begin();
            ConstColIterator colEndIt = M[rowIdx].end();
            ColIterator myColIt = (*this)[rowIdx].begin();
            for (; colIt != colEndIt; ++colIt) {
                while (myColIt.index() < colIt.index())
                    ++ myColIt;
                assert(myColIt.index() == colIt.index());
                
                (*myColIt) = *colIt;
            }
        };
        
        // communicate and add the overlap row contents
        syncronizeEntryValues_(M);
    }
    
private:
    void build_(const BCRSMatrix &M)
    {
        int numDomestic = overlap_->numDomestic();
                       
        // allocate the rows
        this->setSize(numDomestic, numDomestic);
        this->setBuildMode(ParentType::random);

        // communicate the structure of the overlapping rows
        communicateOverlapStructure_(M);

        // copy the rows for the local indices
        int numLocal = M.N();
        for (int rowIdx = 0; rowIdx < numLocal; ++rowIdx)
            this->setrowsize(rowIdx, M.getrowsize(rowIdx));
        
        // set the row size for all additional rows
        typename PeerSet::const_iterator peerIt = overlap_->peerSet().begin();
        typename PeerSet::const_iterator peerEndIt = overlap_->peerSet().end();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            
            typename PeerRows::const_iterator rowIt = overlapEntries_[peerRank].begin();
            typename PeerRows::const_iterator rowEndIt = overlapEntries_[peerRank].end();
            for (; rowIt != rowEndIt; ++rowIt) {              
                int rowIdx = rowIt->first;
                int domesticRowIdx = overlap_->foreignToDomesticIndex(peerRank, rowIdx);
                
                if (overlap_->isFront(peerRank, rowIdx)) {
                    this->setrowsize(domesticRowIdx, 1);
                    continue;
                }

                int nCols = 0;

                // row is a border index, i.e. that we have all the
                // local entries plus the overlap rows from other
                // peers!
                if (overlap_->isLocal(domesticRowIdx))
                    nCols += this->getrowsize(domesticRowIdx);
                
                typename PeerColumns::const_iterator colIt = rowIt->second.begin();
                typename PeerColumns::const_iterator colEndIt = rowIt->second.end();
                for (; colIt != colEndIt; ++colIt) {
                    int colIdx = *colIt;
                    int domesticColIdx = overlap_->foreignToDomesticIndex(peerRank, colIdx);
                
                    // only add new entry if column or row are not on
                    // the border!
                    if (!overlap_->isLocal(domesticColIdx) ||
                        !overlap_->isLocal(domesticRowIdx) )
                        nCols += 1;
                }
                this->setrowsize(domesticRowIdx, nCols);
                assert(this->getrowsize(domesticRowIdx) == nCols);
            }
        }
        this->endrowsizes();

        // add the indices for the local entries
        // copy the rows for the local indices
        for (int rowIdx = 0; rowIdx < numLocal; ++rowIdx) {
            ConstColIterator colIt = M[rowIdx].begin();
            ConstColIterator colEndIt = M[rowIdx].end();
            for (; colIt != colEndIt; ++colIt)
                this->addindex(rowIdx, colIt.index());
        }

        // add the indices for all additional entries
        peerIt = overlap_->peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            
            typename PeerRows::const_iterator rowIt = overlapEntries_[peerRank].begin();
            typename PeerRows::const_iterator rowEndIt = overlapEntries_[peerRank].end();
            for (; rowIt != rowEndIt; ++ rowIt) {
                int rowIdx = rowIt->first;
                int domesticRowIdx = overlap_->foreignToDomesticIndex(peerRank, rowIdx);
               
                if (overlap_->isFront(peerRank, rowIdx)) {
                    this->addindex(domesticRowIdx, domesticRowIdx);
                    continue;
                }                        

                typename PeerColumns::const_iterator colIt = rowIt->second.begin();
                typename PeerColumns::const_iterator colEndIt = rowIt->second.end();
                for (; colIt != colEndIt; ++colIt) {
                    int colIdx = *colIt;
                    int domesticColIdx = overlap_->foreignToDomesticIndex(peerRank, colIdx);
                    
                    // only add new entry if row and column are not a
                    // border index!
                    if (overlap_->isLocal(domesticRowIdx) &&
                        overlap_->isLocal(domesticColIdx))
                    {
                        continue;
                    }

                    this->addindex(domesticRowIdx, domesticColIdx);
                }
            }
        }

        this->endindices();
    };

    // retrieve the number of additional entries for each row in the
    // domestic overlap. this communicates via MPI
    void communicateOverlapStructure_(const BCRSMatrix &M)
    {
        const PeerSet &peerSet = overlap_->peerSet();
        
        int numPeers = peerSet.size();
        int *sendMsgBuff[numPeers][2];
        MPI_Request sendRequests[numPeers][2];

        ///////
        // Send the number of entries from all peers
        ///////
        typename PeerSet::const_iterator it = peerSet.begin();
        typename PeerSet::const_iterator endIt = peerSet.end();
        for (int i = 0; it != endIt; ++it, ++i) {
            int peerRank = *it;

            int msgSize;
            sendMsgBuff[i][0] = createNumEntitiesMsg_(M, msgSize, peerRank);

            // send the size of the overlap of each row to the peer
            // rank
            MPI_Isend(sendMsgBuff[i][0], // pointer to user data 
                      msgSize, // size of user data array
                      MPI_INT, // type of user data
                      peerRank,  // peer rank
                      0, // identifier
                      MPI_COMM_WORLD, // communicator
                      &sendRequests[i][0]); // request object

            // send the overlap indices to the peer rank
            sendMsgBuff[i][1] = createIndicesMsg_(M, msgSize, peerRank);

            MPI_Isend(sendMsgBuff[i][1], // pointer to user data 
                      msgSize, // size of user data array
                      MPI_INT, // type of user data
                      peerRank,  // peer rank
                      0, // identifier
                      MPI_COMM_WORLD, // communicator
                      &sendRequests[i][1]); // request object
        };
        
        ///////
        // Receive all peer indices using syncronous receives
        ///////
        it = peerSet.begin();
        for (; it != endIt; ++it) {
            // unpack the overlap message
            receiveIndices_(*it);
        };

        ///////
        // Wait for all send requests to complete and delete the send
        // buffers
        ///////
        it = peerSet.begin();
        for (int i = 0; it != endIt; ++it, ++i) {
            MPI_Wait(&sendRequests[i][0], MPI_STATUS_IGNORE);
            MPI_Wait(&sendRequests[i][1], MPI_STATUS_IGNORE);

            // free the memory used by the message buffer
            delete[] sendMsgBuff[i][0];
            delete[] sendMsgBuff[i][1];
        }
    }

    // create the message send to a single peer rank by
    // communicateNumEntriesOverlap_()
    int *createNumEntitiesMsg_(const BCRSMatrix &M, int &msgSize, int peerRank)
    {
        msgSize = 
            2*
            overlap_->foreignOverlapWithPeer(peerRank).size();
        int *buff = new int[msgSize];
        const ForeignOverlapWithPeer &olist = overlap_->foreignOverlapWithPeer(peerRank);
        typename ForeignOverlapWithPeer::const_iterator it = olist.begin();
        typename ForeignOverlapWithPeer::const_iterator endIt = olist.end();
        for (int i = 0; it != endIt; ++it, ++i) {
            int rowIdx = std::get<0>(*it);
            buff[2*i + 0] = rowIdx;
            if (overlap_->isForeignFront(peerRank, rowIdx))
                buff[2*i + 1] = 0; // we do not send anything for front rows
            else
                buff[2*i + 1] = M.getrowsize(rowIdx);
        }
        return buff;
    };

    // create the message send to a single peer rank by
    // communicateNumEntriesOverlap_()
    int *createIndicesMsg_(const BCRSMatrix &M, int &msgSize, int peerRank)
    {
        msgSize = 0;

        const ForeignOverlapWithPeer &olist = overlap_->foreignOverlapWithPeer(peerRank);
        
        // calculate message size
        typename ForeignOverlapWithPeer::const_iterator it = olist.begin();
        typename ForeignOverlapWithPeer::const_iterator endIt = olist.end();
        for (; it != endIt; ++it) {
            int rowIdx = std::get<0>(*it);
            msgSize += M.getrowsize(rowIdx);
        }

        // allocate message buffer
        int *buff = new int[msgSize];
        msgSize = 0;

        // fill the message buffer
        it = olist.begin();
        for (int i = 0; it != endIt; ++it) {
            int rowIdx = std::get<0>(*it);
            
            if (overlap_->isForeignFront(peerRank, rowIdx))
                continue; // we do not send anything for front rows

            ConstColIterator colIt = M[rowIdx].begin();
            ConstColIterator colEndIt = M[rowIdx].end();
            for (; colIt != colEndIt; ++colIt, ++i) {
                buff[i] = colIt.index();
                ++ msgSize;
            }
        }
        
        return buff;
    };

    // retrieve and process the message by a single peer rank. called
    // from communicateNumEntriesOverlap_()
    void receiveIndices_(int peerRank)
    {
        int numRows = overlap_->domesticOverlapWithPeer(peerRank).size();
        overlapEntries_[peerRank].resize(numRows);

        numOverlapEntries_[peerRank] = 0;
        int msgSize = 2*numRows;
        int *buff = new int[msgSize];
        // receive size of overlap message buffer
        MPI_Recv(buff, // receive message buffer
                 msgSize, // receive message size
                 MPI_INT, // object type
                 peerRank, // peer rank
                 0, // identifier
                 MPI_COMM_WORLD, // communicator
                 MPI_STATUS_IGNORE); // status

        for (int i = 0; i < numRows; ++i) {
            int rowIdx = buff[2*i + 0];
            int numCols = buff[2*i + 1];

            overlapEntries_[peerRank][i].first = rowIdx;
            overlapEntries_[peerRank][i].second.resize(numCols);
            numOverlapEntries_[peerRank] += numCols;
        }
        delete[] buff;
        
        // receive overlap indices themselfs
        msgSize = numOverlapEntries_[peerRank];
        buff = new int[msgSize];
        MPI_Recv(buff, // receive message buffer
                 msgSize, // receive message size
                 MPI_INT, // object type
                 peerRank, // peer rank
                 0, // identifier
                 MPI_COMM_WORLD, // communicator
                 MPI_STATUS_IGNORE); // status

        int pos = 0;
        for (int i = 0; i < numRows; ++i) {
            int rowIdx = overlapEntries_[peerRank][i].first;
            int nCols = overlapEntries_[peerRank][i].second.size();
            for (int j = 0; j < nCols; ++j) {
                int colIndex = buff[pos];
                ++pos;
                overlapEntries_[peerRank][i].second[j] = colIndex;
            }

            if (overlap_->isFront(peerRank, rowIdx))
                overlapEntries_[peerRank][i].second.resize(0);
        }
        delete[] buff;
    };
    
    void syncronizeEntryValues_(const BCRSMatrix &M)
    {
        const PeerSet &peerSet = overlap_->peerSet();
        
        int numPeers = peerSet.size();
        field_type *sendMsgBuff[numPeers];
        MPI_Request sendRequests[numPeers];

        ///////
        // Send the number of entries from all peers
        ///////
        typename PeerSet::const_iterator it = peerSet.begin();
        typename PeerSet::const_iterator endIt = peerSet.end();
        for (int i = 0; it != endIt; ++it, ++i) {
            int peerRank = *it;

            // send the overlap indices to the peer rank
            int msgSize;
            sendMsgBuff[i] = createValuesMsg_(M, msgSize, peerRank);

#warning MPI_DOUBLE is only correct if field_type == double!
            MPI_Isend(sendMsgBuff[i], // pointer to user data 
                      msgSize, // size of user data array
                      MPI_DOUBLE, // type of user data
                      peerRank,  // peer rank
                      0, // identifier
                      MPI_COMM_WORLD, // communicator
                      &sendRequests[i]); // request object
        };
        
        ///////
        // Receive all peer indices using syncronous receives
        ///////
        it = peerSet.begin();
        for (; it != endIt; ++it) {
            // unpack the overlap message
            receiveValues_(*it);
        };

        ///////
        // Wait for all send requests to complete and delete the send
        // buffers
        ///////
        it = peerSet.begin();
        for (int i = 0; it != endIt; ++it, ++i) {
            MPI_Wait(&sendRequests[i], MPI_STATUS_IGNORE);

            // free the memory used by the message buffer
            delete[] sendMsgBuff[i];
        }
    };

    field_type *createValuesMsg_(const BCRSMatrix &M, int &msgSize, int peerRank)
    {
        msgSize = 0;

        const ForeignOverlapWithPeer &olist = overlap_->foreignOverlapWithPeer(peerRank);
        
        // calculate message size
        typename ForeignOverlapWithPeer::const_iterator it = olist.begin();
        typename ForeignOverlapWithPeer::const_iterator endIt = olist.end();
        for (; it != endIt; ++it) {
            int rowIdx = std::get<0>(*it);
            msgSize += M.getrowsize(rowIdx);
        }
        msgSize *= block_type::rows*block_type::cols;

        // allocate message buffer
        field_type *buff = new field_type[msgSize];
        int pos = 0;

        // fill the message buffer
        it = olist.begin();
        for (int i = 0; it != endIt; ++it) {
            int rowIdx = std::get<0>(*it);
            
            if (overlap_->isForeignFront(peerRank, rowIdx))
                continue; // we do not send anything for front rows
            
            ConstColIterator colIt = M[rowIdx].begin();
            ConstColIterator colEndIt = M[rowIdx].end();
            for (; colIt != colEndIt; ++colIt, ++i) {
                for (int i = 0; i < block_type::rows; ++i) {
                    for (int j = 0; j < block_type::cols; ++j) {
                        buff[pos] = (*colIt)[i][j];
                        ++pos;
                    }
                }
            }
        }
        assert(pos <= msgSize);
        msgSize = pos;

        return buff;
    };

    void receiveValues_(int peerRank)
    {
        // receive the values of the overlap entries
        // calculate message size
        const DomesticOverlapWithPeer &olist = overlap_->domesticOverlapWithPeer(peerRank);
        int numRows = olist.size();
        int msgSize = numOverlapEntries_[peerRank] * block_type::rows*block_type::cols;
        field_type *buff = new field_type[msgSize];

#warning MPI_DOUBLE is only correct if field_type == double!
        MPI_Recv(buff, // receive message buffer
                 msgSize, // receive message size
                 MPI_DOUBLE, // object type
                 peerRank, // peer rank
                 0, // identifier
                 MPI_COMM_WORLD, // communicator
                 MPI_STATUS_IGNORE); // status

        int pos = 0;
        for (int i = 0; i < numRows; ++i) {
            int rowIdx = overlapEntries_[peerRank][i].first;
            int domRowIdx = overlap_->foreignToDomesticIndex(peerRank, rowIdx);

            if (overlap_->isFront(peerRank, rowIdx)) {
                // the front entries only get an identity matrix
                block_type &x = (*this)[domRowIdx][domRowIdx];
                for (int k = 0; k < block_type::rows; ++k)
                    x[k][k] = 1.0;
                // front rows do not send anyting
                continue;
            }

            int nCols = overlapEntries_[peerRank][i].second.size();
            for (int j = 0; j < nCols; ++j) {
                int colIdx = overlapEntries_[peerRank][i].second[j];
                int domColIdx = overlap_->foreignToDomesticIndex(peerRank, colIdx);
                     
                block_type &x = (*this)[domRowIdx][domColIdx];
                for (int k1 = 0; k1 < block_type::rows; ++k1) {
                    for (int k2 = 0; k2 < block_type::cols; ++k2) {
                        x[k1][k2] += buff[pos];
                        assert(std::isfinite(buff[pos]));
                        ++pos;
                    }
                }
            }
        }
        assert(pos == msgSize);
        delete[] buff;
    };

    const Overlap *overlap_;
    std::map<ProcessRank, int> numOverlapEntries_;
    OverlapStructure overlapEntries_;
};

template <class BlockVector, class Overlap>
class OverlapBlockVector : public BlockVector
{
    typedef BlockVector ParentType;
    typedef OverlapBlockVector<BlockVector, Overlap> ThisType;

    typedef typename Overlap::PeerSet PeerSet;
    typedef typename Overlap::ForeignOverlapWithPeer ForeignOverlapWithPeer;
    typedef typename Overlap::DomesticOverlapWithPeer DomesticOverlapWithPeer;
    

public:
    typedef typename ParentType::field_type field_type;
    typedef typename ParentType::block_type block_type;
    
    OverlapBlockVector(const BlockVector &bv, const Overlap &overlap)
        : ParentType(overlap.numDomestic()),
          overlap_(&overlap)
    {
        set(bv);
    }

    // copy constructor
    OverlapBlockVector(const ThisType &v)
        : ParentType(v), overlap_(v.overlap_)
    { 
    };

    void set(const BlockVector &v)
    {
        // copy local entries
        int numLocal = overlap_->numLocal();
        for (int i = 0; i < numLocal; ++i)
            (*this)[i] = v[i];

        // set the additional DOFs to 0
        int numDomestic = overlap_->numDomestic();
        for (int i = numLocal; i < numDomestic; ++i)
            (*this)[i] = 0.0;
        
        // retrieve the overlap values from their respective owners
        // copyOverlapFromOwner();
    }
    
    void markAdditional()
    {
        // set the additional DOFs to 0
        int numLocal = overlap_->numLocal();
        int numDomestic = overlap_->numDomestic();
        for (int i = numLocal; i < numDomestic; ++i)
            (*this)[i] = 123.0;
    };

    void assignToNonOverlapping(BlockVector &v)
    {
        int numLocal = overlap_->numLocal();
        for (int i = 0; i < numLocal; ++i)
            v[i] = (*this)[i];
    }

    void weightedSumOnOverlap()
    {
        const PeerSet &peerSet = overlap_->peerSet();
        
        int numPeers = peerSet.size();
        field_type *sendMsgBuff[numPeers];
        MPI_Request sendRequests[numPeers];

        ///////
        // Send the values of all own entries to all peers
        ///////
        typename PeerSet::const_iterator it = peerSet.begin();
        typename PeerSet::const_iterator endIt = peerSet.end();
        for (int i = 0; it != endIt; ++it, ++i) {
            int peerRank = *it;

            int msgSize;
            sendMsgBuff[i] = createValuesMsg_(msgSize, peerRank);
            
            // send the size of the overlap of each row to the peer
            // rank
#warning MPI_DOUBLE is only correct if field_type == double!
            MPI_Isend(sendMsgBuff[i], // pointer to user data 
                      msgSize, // size of user data array
                      MPI_DOUBLE, // type of user data
                      peerRank,  // peer rank
                      0, // identifier
                      MPI_COMM_WORLD, // communicator
                      &sendRequests[i]); // request object

        };

        // divide all values by the number of processes which "see" a
        // given index
        int numDomestic = overlap_->numDomestic();
        for (int i = 0; i < numDomestic; ++i) {
            // the number of processes which "see" the row is the
            // number of peer processes where the row is in the
            // overlap plus ourselfs
            int nProc = 1 + overlap_->numOverlapPeers(i);
            (*this)[i] /= nProc;
        };

        ///////
        // Receive all peer values using syncronous receives
        ///////
        it = peerSet.begin();
        for (; it != endIt; ++it) {
            // unpack the overlap message
            receiveValuesWeightedSum_(*it);
        };

        ///////
        // Wait for all send requests to complete and delete the send
        // buffers
        ///////
        it = peerSet.begin();
        for (int i = 0; it != endIt; ++it, ++i) {
            MPI_Wait(&sendRequests[i], MPI_STATUS_IGNORE);

            // free the memory used by the message buffer
            delete[] sendMsgBuff[i];
        }
        
    };
    
    void copyOverlapFromOwner()
    {
        const PeerSet &peerSet = overlap_->peerSet();
        
        int numPeers = peerSet.size();
        field_type *sendMsgBuff[numPeers];
        MPI_Request sendRequests[numPeers];

        ///////
        // Send the values of all own entries to all peers
        ///////
        typename PeerSet::const_iterator it = peerSet.begin();
        typename PeerSet::const_iterator endIt = peerSet.end();
        for (int i = 0; it != endIt; ++it, ++i) {
            int peerRank = *it;

            int msgSize;
            sendMsgBuff[i] = createValuesMsg_(msgSize, peerRank);
            
            // send the size of the overlap of each row to the peer
            // rank
#warning MPI_DOUBLE is only correct if field_type == double!
            MPI_Isend(sendMsgBuff[i], // pointer to user data 
                      msgSize, // size of user data array
                      MPI_DOUBLE, // type of user data
                      peerRank,  // peer rank
                      0, // identifier
                      MPI_COMM_WORLD, // communicator
                      &sendRequests[i]); // request object

        };

        ///////
        // Receive all peer indices using syncronous receives
        ///////
        it = peerSet.begin();
        for (; it != endIt; ++it) {
            // unpack the overlap message
            receiveValues_(*it);
        };

        ///////
        // Wait for all send requests to complete and delete the send
        // buffers
        ///////
        it = peerSet.begin();
        for (int i = 0; it != endIt; ++it, ++i) {
            MPI_Wait(&sendRequests[i], MPI_STATUS_IGNORE);

            // free the memory used by the message buffer
            delete[] sendMsgBuff[i];
        }
    };

    ThisType &operator=(const ThisType &v)
    {
        ParentType::operator=(v);
        return *this;
    }

    ThisType &operator=(const field_type v)
    {
        ParentType::operator=(v);
        return *this;
    }
    
private:
    field_type *createValuesMsg_(int &msgSize, int peerRank)
    {
        const ForeignOverlapWithPeer &olist = overlap_->foreignOverlapWithPeer(peerRank);
        int overlapSize = olist.size();

        // calculate message size
        msgSize = block_type::size * overlapSize;

        // allocate message buffer
        field_type *buff = new field_type[msgSize];
        int pos = 0;

        // fill the message buffer
        typename ForeignOverlapWithPeer::const_iterator it = olist.begin();
        typename ForeignOverlapWithPeer::const_iterator endIt = olist.end();
        for (; it != endIt; ++it) {
            int rowIdx = std::get<0>(*it);
            
            for (int i = 0; i < block_type::size; ++i) {
                buff[pos] = (*this)[rowIdx][i];
                ++pos;
            }
        }

        return buff;
    };

    void receiveValues_(int peerRank)
    {
        // receive overlap indices itself
        const DomesticOverlapWithPeer &olist = overlap_->domesticOverlapWithPeer(peerRank);
        int overlapSize = olist.size();
        int msgSize = block_type::size * overlapSize;
        field_type *buff = new field_type[msgSize];

#warning MPI_DOUBLE is only correct if field_type == double!
        MPI_Recv(buff, // receive message buffer
                 msgSize, // receive message size
                 MPI_DOUBLE, // object type
                 peerRank, // peer rank
                 0, // identifier
                 MPI_COMM_WORLD, // communicator
                 MPI_STATUS_IGNORE); // status

        int pos = 0;
        typename DomesticOverlapWithPeer::const_iterator it = olist.begin();
        typename DomesticOverlapWithPeer::const_iterator endIt = olist.end();
        for (; it != endIt; ++it) {
            int rowIdx = std::get<0>(*it);
            int domRowIdx = overlap_->foreignToDomesticIndex(peerRank, rowIdx);
            for (int i = 0; i < block_type::size; ++i) {
                (*this)[domRowIdx][i] = buff[pos];
                ++pos;
            }

            /*if (overlap_->isLocal(domRowIdx)) {
                std::cout << "border @" << domRowIdx << "\n";
                (*this)[domRowIdx] = 123;
            }
            */
        }

        delete[] buff;

    };

    void receiveValuesWeightedSum_(int peerRank)
    {
        // receive overlap indices itself
        const DomesticOverlapWithPeer &olist = overlap_->domesticOverlapWithPeer(peerRank);
        int overlapSize = olist.size();
        int msgSize = block_type::size * overlapSize;
        field_type *buff = new field_type[msgSize];

#warning MPI_DOUBLE is only correct if field_type == double!
        MPI_Recv(buff, // receive message buffer
                 msgSize, // receive message size
                 MPI_DOUBLE, // object type
                 peerRank, // peer rank
                 0, // identifier
                 MPI_COMM_WORLD, // communicator
                 MPI_STATUS_IGNORE); // status

        int pos = 0;

        typename DomesticOverlapWithPeer::const_iterator it = olist.begin();
        typename DomesticOverlapWithPeer::const_iterator endIt = olist.end();
        for (; it != endIt; ++it) {
            int rowIdx = std::get<0>(*it);
            int domRowIdx = overlap_->foreignToDomesticIndex(peerRank, rowIdx);
            
            // the number of processes which "see" the row is the
            // number of processes where the row is in the overlap
            // plus ourselfs
            int numProc = 1 + overlap_->numOverlapPeers(domRowIdx);
            
            for (int i = 0; i < block_type::size; ++i) {
                (*this)[domRowIdx][i] += buff[pos] / numProc;
                ++pos;
            }
        }

        delete[] buff;

    };

    const Overlap *overlap_;
};

template <class BlockVector, class Overlap>
class OverlapScalarProduct : public Dune::ScalarProduct<BlockVector>
{
public:
    typedef typename BlockVector::field_type field_type;
    typedef BlockVector domain_type;

#warning should probably be overlapping!?
    enum { category = Dune::SolverCategory::sequential };
    
    OverlapScalarProduct(const Overlap &overlap)
        : overlap_(&overlap)
    {};
    
    field_type dot(const BlockVector &x, const BlockVector &y)
    { 
        field_type localRes = 0;
        int n = overlap_->numLocal();
        for (int i = 0; i < n; ++i) {
            localRes += x[i]*y[i];
        };
        
        field_type globalRes = 0;
#warning MPI_DOUBLE is only correct if field_type == double!
        MPI_Allreduce(&localRes, // source buffer
                      &globalRes, // destination buffer
                      1, // number of objects in buffers
                      MPI_DOUBLE, // data type
                      MPI_SUM, // operation
                      MPI_COMM_WORLD); // communicator
        return globalRes;
    };

    field_type norm(const BlockVector &x)
    { 
        field_type localRes = 0;
        int n = overlap_->numLocal();
        for (int i = 0; i < n; ++i) {
            localRes += x[i]*x[i];
        };
        
        field_type globalRes = 0;
#warning MPI_DOUBLE is only correct if field_type == double!
        MPI_Allreduce(&localRes, // source buffer
                      &globalRes, // destination buffer
                      1, // number of objects in buffers
                      MPI_DOUBLE, // data type
                      MPI_SUM, // operation
                      MPI_COMM_WORLD); // communicator
        return std::sqrt(globalRes);
    };

private:
    const Overlap *overlap_;
};

template <class SeqPreCond, class OverlapMatrix, class Overlap>
class OverlapPreconditioner : 
    public Dune::Preconditioner<typename SeqPreCond::domain_type, 
                                typename SeqPreCond::range_type>
{
public:
    typedef typename SeqPreCond::domain_type domain_type;
    typedef typename SeqPreCond::range_type range_type;
    typedef typename SeqPreCond::field_type field_type;

#warning should probably be overlapping!?
    enum { category = Dune::SolverCategory::sequential };

    OverlapPreconditioner(SeqPreCond &seqPreCond,
                          const Overlap &overlap)
        : seqPreCond_(seqPreCond), overlap_(&overlap)
    {
    }

    void pre(domain_type &x, range_type &y)
    {
        seqPreCond_.pre(x, y);
    };

    void apply(domain_type &x, const range_type &y)
    {
        seqPreCond_.apply(x, y);

        // communicate the results on the overlap
        x.weightedSumOnOverlap();
    };
   
    void post(domain_type &x)
    {
        seqPreCond_.post(x);
    };

private:
    SeqPreCond seqPreCond_;
    const Overlap *overlap_;
};

} // namespace Dumux

#endif
