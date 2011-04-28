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
    typedef std::map<ProcessRank, std::list<IndexDistanceNpeers> > ForeignOverlapByRank;
    typedef std::vector<std::map<ProcessRank, BorderDistance> > ForeignOverlapByIndex;
    typedef std::map<ProcessRank, std::vector<IndexDistanceNpeers> > DomesticOverlap;
    typedef std::map<ProcessRank, std::set<PeerIndex> > DomesticFront;
    typedef std::map<ProcessRank, std::map<PeerIndex, DomesticIndex> > DomesticBorder;

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
    }

    bool isLocal(int idx) const
    { return idx < numLocalIndices_; };

    bool isBorder(int peerRank, int foreignIdx) const
    { return domesticBorder_.find(peerRank)->second.count(foreignIdx) > 0; };

    bool isFront(int peerRank, int foreignIdx) const
    { return domesticFront_[peerRank].count(foreignIdx) > 0; };

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

        if (isLocal(domesticIdx)) {
            // we do not map the border indices back to foreign
            // indices!
            return -1;
        }

        int offset =
            domesticIdx -
            numLocalIndices_ -
            domesticOffset_.find(peerRank)->second;
        const IndexDistanceNpeers &tmp = 
            domesticOverlap_.find(peerRank)->second[offset];
        return std::get<0>(tmp);
    };

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
        for (; it != endIt; ++it) {
            int peerRank = *it;
            assert(peerRank != myRank);

            int sendMsgSize;
            sendMsgBuff[peerRank] = packSendMessage_(peerRank, sendMsgSize);
            
            //std::cout << "rank " << myRank << " sends overlap message of size " << sendMsgSize << " to rank " << *it << "\n";

            // send the size of the overlap message buffer to the
            // peer rank
            MPI_Isend(&sendMsgSize, // pointer to user data 
                      1, // size of user data array
                      MPI_INT, // type of user data
                      peerRank,  // peer rank
                      0, // identifier
                      MPI_COMM_WORLD, // communicator
                      &sendRequests[*it][0]); // request object
            // send the overlap indices and distances to the peer
            // rank
            MPI_Isend(sendMsgBuff[peerRank], // pointer to user data 
                      sendMsgSize, // size of user data array
                      MPI_INT, // type of user data
                      peerRank,  // peer rank
                      0, // identifier
                      MPI_COMM_WORLD, // communicator
                      &sendRequests[*it][1]); // request object
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
            
            //std::cout << "rank " << myRank << " received overlap message of size " << recvMsgSize << " from rank " << *it << "\n";
        };

        ///////
        // Wait for all send requests to complete, (this should not be
        // necessary, since we used syncronous receives, but it might
        // avoid memory leaks inside the MPI library) and delete the
        // send buffers
        ///////
        it = peerSet_.begin();
        for (; it != endIt; ++it) {
            MPI_Wait(&sendRequests[*it][0], MPI_STATUS_IGNORE);
            MPI_Wait(&sendRequests[*it][1], MPI_STATUS_IGNORE);

            // free the memory used by the message buffer
            delete[] sendMsgBuff[*it];
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
        };

        return buff;
    };

    // unpack the overlap message buffer received from a peer
    void unpackRecvMessage_(int peerRank, int *overlapMsg, int &msgSize) 
    {
        int n = msgSize / 3;
        
        // create the domestic overlap for the peer rank
        for (int i = 0; i < n; ++i) {
            int peerIndexIdx = overlapMsg[3*i + 0];
            int peerDistance = overlapMsg[3*i + 1];
            int numRanks = overlapMsg[3*i + 2];
            assert(peerRank >= 0);
            assert(peerIndexIdx >= 0);
            IndexDistanceNpeers tmp(peerIndexIdx, peerDistance, numRanks);
            domesticOverlap_[peerRank].push_back(tmp);
        };
    };

    // calculate the border indices given the initial seed list
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
        };
    };
    
    // builds the map to translate from a (peer rank, foreign index)
    // pair to the domestic index
    void buildForeignToDomesticMap_()
    {
        // calculate the offsets for the additional domestic indices
        int n = 0;
        DomesticOverlap::iterator ovrlpIt = domesticOverlap_.begin();
        DomesticOverlap::iterator ovrlpEndIt = domesticOverlap_.end();
        for (; ovrlpIt != ovrlpEndIt; ++ovrlpIt) {
            int peerRank = ovrlpIt->first;

            // calculate the number of border indices for the current
            // peer rank.
            int numBorder = 0;
            int numIndices = ovrlpIt->second.size();
            for (int i = 0; i < numIndices; ++i) {
                int borderDist = std::get<1>(ovrlpIt->second[i]);

                if (borderDist == 0) {
                    // move this index to the back
                    std::swap(ovrlpIt->second[i],
                              ovrlpIt->second[numIndices - 1 - numBorder]);

                    // increment border index count
                    ++ numBorder;
                    continue;
                }
            }
                        
            //borderIndices_[peerRank] = nBorder;
            domesticOffset_[peerRank] = n;
            n += numIndices - numBorder;
        }

        // add the border rows to the map
        DomesticBorder::const_iterator brdrIt = domesticBorder_.begin();
        DomesticBorder::const_iterator brdrEndIt = domesticBorder_.end();
        for (; brdrIt != brdrEndIt; ++brdrIt) {
            int peerRank = brdrIt->first;
            std::map<PeerIndex, DomesticIndex>::const_iterator it = brdrIt->second.begin();
            std::map<PeerIndex, DomesticIndex>::const_iterator endIt = brdrIt->second.end();
            for (; it != endIt; ++it) {
                int localIdx = it->first;
                int peerIdx = it->second;
                foreignToDomesticMap_[peerRank][peerIdx] = localIdx;
            }
        };        

        // add the additional indices to the {(foreignIndex, peerRank)
        // -> domestic index} map
        ovrlpIt = domesticOverlap_.begin();
        for (; ovrlpIt != ovrlpEndIt; ++ovrlpIt) {
            int peerRank = ovrlpIt->first;

            int numIndices = ovrlpIt->second.size();
            for (int i = 0; i < numIndices; ++i) {
                int peerIdx = std::get<0>(ovrlpIt->second[i]);
                int borderDist = std::get<1>(ovrlpIt->second[i]);
                int domIdx = numLocalIndices_ + domesticOffset_[peerRank] + i;

                if (borderDist == 0) {
                    // ignore border indices since they are not
                    // additional degrees of freedom
                    continue;
                }

                foreignToDomesticMap_[peerRank][peerIdx] = domIdx;
            };
        };
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

    // stores the "border" indices of the domestic overlap
    DomesticBorder domesticBorder_;

    // number of local indices
    int numLocalIndices_;
    
    // offset for the overlap indices for each peer process in number
    // of indices after the last "owned" index
    std::map<ProcessRank, int> domesticOffset_;

    // maps a peer local index to a domestic local index
    std::map<ProcessRank, std::map<PeerIndex, DomesticIndex> > foreignToDomesticMap_;
};

} // namespace Dumux

#endif
