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
class VertexSeedListFromGrid : public Dune::CommDataHandleIF< VertexSeedListFromGrid<GridView, VertexMapper>,
                                                              int >
{
    typedef int ProcessRank;
    typedef int RowIndex;
    typedef std::pair<RowIndex, ProcessRank> RowProcessPair;
    typedef std::list<RowProcessPair> SeedList;

public:
    VertexSeedListFromGrid(const GridView &gv, 
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
    { return 1; };

    template<class MessageBufferImp, class EntityType> 
    void gather(MessageBufferImp &buff, const EntityType &e) const 
    { buff.write(gv_.comm().rank()); };

    template<class MessageBufferImp, class EntityType> 
    void scatter(MessageBufferImp &buff, const EntityType &e, size_t n)
    {
        int remoteRank;
        buff.read(remoteRank);
        int vertIdx = map_.map(e);
        seedList_.push_back(RowProcessPair(vertIdx, remoteRank));
    };
    
    // Access to the initial seed list.
    const SeedList &seedList() const
    { return seedList_; }

private:
    const GridView gv_;
    const VertexMapper &map_;
    SeedList seedList_;
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
    typedef int RowIndex;
    typedef std::pair<RowIndex, ProcessRank> RowProcessPair;
    typedef std::pair<RowIndex, BorderDistance> RowDistancePair;
    typedef std::list<RowProcessPair> SeedList;
    typedef std::pair<ProcessRank, int> ProcessDistancePair;

    typedef std::set<ProcessRank> PeerSet;
    typedef std::map<ProcessRank, std::list<RowDistancePair> > ForeignOverlapByRank;
    typedef std::vector<std::map<ProcessRank, BorderDistance> > ForeignOverlapByRow;
    typedef std::map<ProcessRank, std::vector<RowDistancePair> > DomesticOverlap;
    typedef std::map<ProcessRank, std::list<RowIndex> > DomesticFront;
    typedef std::map<ProcessRank, std::list<RowIndex> > DomesticBorder;

    OverlapFromBCRSMatrix(const BCRSMatrix &A,
                          const SeedList &seedList,
                          int overlapSize)
    {
        // find the set of processes which have an overlap with the
        // local processes. (i.e. the set of processes which we will
        // have to communicate to.)
        buildPeerSet_(seedList);

        // calculate the foreign overlap for the local partition,
        // i.e. find the distance of each row from the seed set.
        foreignOverlapByRow_.resize(A.N());
        extendForeignOverlap_(A, seedList, 0, overlapSize);

        // group foreign overlap by remote process rank
        groupForeignOverlapByRank_(); 
        
        // calculate the domestic overlap (i.e. all overlap indices in
        // foreign processes which the current process overlaps.)
        // This requires communication via MPI.
        buildDomesticOverlap_();

        // calculate the border and front lists
        buildBorderFront_(overlapSize);

        // set the number of local indices to the number of rows of
        // the BCRS matrix
        numDomesticIndices_ = A.N();

        // calculate the per-process offsets for the additional
        // indices due to the overlap.
        calculateOverlapOffsets_();
    }

    bool isLocal(int idx) const
    { return idx < numDomesticIndices_; };

    int foreignToDomesticIndex(int peerRank, int foreignIdx) const
    { 
        assert(foreignToLocalMap_.count(peerRank) > 0);
        assert(foreignToLocalMap_[peerRank].count(foreignIdx) > 0);

        return foreignToLocalMap_[peerRank][foreignIdx];
    };

    int domesticToForeignIndex(int domesticIdx, int peerRank) const
    { 
        assert(domesticOffset_.count(peerRank) > 0);
        assert(domesticIdx >= numDomesticIndices_ + domesticOffset_[peerRank]);
        assert(domesticIdx < numDomesticIndices_ + domesticOffset_[peerRank] + domesticOverlap_[peerRank].size());

        return domesticOverlap_[peerRank][domesticIdx - numDomesticIndices_ - domesticOffset_[peerRank]];
    };

    void printOverlap() const
    {
        int myRank;
        MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
        std::cout << "Number of ranks with overlap: " << peerSet_.size() << "\n";

        DomesticOverlap::const_iterator it = domesticOverlap_.begin();
        DomesticOverlap::const_iterator endIt = domesticOverlap_.end();
        for (; it != endIt; ++it) {
            std::cout << myRank << " overlaps with rank " << it->first << " at rows(distance): " << ": ";
            
            std::vector<RowDistancePair>::const_iterator rowIt = it->second.begin();
            std::vector<RowDistancePair>::const_iterator rowEndIt = it->second.end();
            for (; rowIt != rowEndIt; ++rowIt) {
                std::cout << rowIt->first << "(" << rowIt->second << ") ";
            };
            std::cout << "\n";
        }

/*
        ForeignOverlapByRank::const_iterator it = foreignOverlapByRank_.begin();
        ForeignOverlapByRank::const_iterator endIt = foreignOverlapByRank_.end();
        for (; it != endIt; ++it) {
            std::cout << "Overlap rows(distance) for rank " << it->first << ": ";
            
            std::list<RowDistancePair>::const_iterator rowIt = it->second.begin();
            std::list<RowDistancePair>::const_iterator rowEndIt = it->second.end();
            for (; rowIt != rowEndIt; ++rowIt) {
                std::cout << rowIt->first << "(" << rowIt->second << ") ";
            };
            std::cout << "\n";
        }
*/
/*
        int n = foreignOverlapByRow_.size();
        for (int i = 0; i < n; ++i) {
            const std::map<int, int> &rowSet = foreignOverlapByRow_[i];
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
                               const SeedList &seedList,
                               int overlapDistance,
                               int overlapSize)
    {
        // add all processes in the seed rows of the current overlap
        // level
        SeedList::const_iterator it = seedList.begin();
        SeedList::const_iterator endIt = seedList.end();
        for (; it != endIt; ++it) {
            if (foreignOverlapByRow_[it->first].count(it->second) == 0)
                foreignOverlapByRow_[it->first][it->second] = overlapDistance;
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
                if (foreignOverlapByRow_[newIdx].count(processIdx) > 0)
                    continue;
                    
                // check whether the new index is already in the next seed list
                RowProcessPair newPair(newIdx, processIdx);
                if (std::find(nextSeedList.begin(), 
                              nextSeedList.end(),
                              newPair) != nextSeedList.end())
                    continue; // we already have this pair                    

                // add the current processes to the seed list for the
                // next overlap level
                nextSeedList.push_back(newPair);
            }
        }
            
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
        int nRows = foreignOverlapByRow_.size();
        for (int i = 0; i < nRows; ++i)
        {
            std::map<ProcessRank, BorderDistance>::const_iterator it = 
                foreignOverlapByRow_[i].begin();
            std::map<ProcessRank, BorderDistance>::const_iterator endIt = 
                foreignOverlapByRow_[i].end();
            for (; it != endIt; ++it) 
                foreignOverlapByRank_[it->first].push_back(RowDistancePair(i, 
                                                                           it->second));
        };
    }
    
    // build the domestic overlap of the process, i.e. the overlapping
    // indices of the current process within the peer processes
    void buildDomesticOverlap_()
    {
        // communicate the overlapping indices to the remote processes:
        
        
        ///////
        // Send all indices using asyncronous send
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
            // remote rank
            MPI_Isend(&sendMsgSize, // pointer to user data 
                      1, // size of user data array
                      MPI_INT, // type of user data
                      peerRank,  // remote rank
                      0, // identifier
                      MPI_COMM_WORLD, // communicator
                      &sendRequests[*it][0]); // request object
            // send the overlap indices and distances to the remote
            // rank
            MPI_Isend(sendMsgBuff[peerRank], // pointer to user data 
                      sendMsgSize, // size of user data array
                      MPI_INT, // type of user data
                      peerRank,  // remote rank
                      0, // identifier
                      MPI_COMM_WORLD, // communicator
                      &sendRequests[*it][1]); // request object
        };

        ///////
        // Receive all remote indices using syncronous receives
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
        // Wait for all send requests to complete. (this should not be
        // necessary, since we used syncronous receives, but it might
        // avoid memory leaks inside the MPI library.)
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
        msgSize = 2*n;
        int *buff = new int[msgSize];

        std::list<RowDistancePair>::const_iterator it = foreignOverlapByRank_[peerRank].begin();
        std::list<RowDistancePair>::const_iterator endIt = foreignOverlapByRank_[peerRank].end();
        for (int i = 0; i < n; ++i, ++it) {
            assert(it->first >= 0);
            assert(it->second >= 0);
            buff[2*i + 0] = it->first; // row index
            buff[2*i + 1] = it->second; // border distance
        };

        return buff;
    };

    // unpack the overlap message buffer received from a peer
    void unpackRecvMessage_(int peerRank, int *overlapMsg, int &msgSize) 
    {
        int n = msgSize / 2;
        
        // create the domestic overlap for the peer rank
        for (int i = 0; i < n; ++i) {
            int peerRowIdx = overlapMsg[2*i + 0];
            int peerDistance = overlapMsg[2*i + 1];
            assert(peerRank >= 0);
            assert(peerRowIdx >= 0);
            domesticOverlap_[peerRank].push_back(RowDistancePair(peerRowIdx, peerDistance));
        };
    };

    // assuming that the domestic overlap has already been calulated,
    // calculate the "border" and "front" indices from the overlap
    // distance
    void buildBorderFront_(int overlapSize)
    {
        int myRank;
        MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

        std::set<int>::iterator it = peerSet_.begin();
        std::set<int>::iterator endIt = peerSet_.end();
        for (; it != endIt; ++it) {
            int peerRank = *it;
            assert(peerRank != myRank);
            
            std::vector<RowDistancePair >::const_iterator oit = domesticOverlap_[peerRank].begin();
            std::vector<RowDistancePair >::const_iterator oEndIt = domesticOverlap_[peerRank].end();
            for (; oit != oEndIt; ++oit) {
                int peerRowIdx = oit->first;
                int peerDistance = oit->second;
                if (peerDistance == 0) 
                    // border index
                    domesticOverlapBorder_[peerRank].push_back(peerRowIdx);
                if (peerDistance == overlapSize) 
                    // front index. TODO: this is only true for
                    // indices which are not on the boundary of the
                    // remote proces' domain. we probably have to
                    // communicate to get the real front set!
                    domesticOverlapBorder_[peerRank].push_back(peerRowIdx);
            }
        };
    };

    // calculate the offsets for the additional indices
    void calculateOverlapOffsets_()
    {
        std::set<int>::iterator it = peerSet_.begin();
        std::set<int>::iterator endIt = peerSet_.end();
        int n = 0;
        for (; it != endIt; ++it) {
            int peerRank = *it;
            domesticOffset_[peerRank] = n;
            n += domesticOverlap_[peerRank].size();
        }
    };

    // set of processes with which we have to communicate
    PeerSet peerSet_;
    
    // stores the set of process ranks which are in the overlap for a
    // given row index "owned" by the current rank. The second value
    // store the distance from the nearest process border.
    ForeignOverlapByRow foreignOverlapByRow_;

    // stores a list of foreign overlap indices for each rank
    ForeignOverlapByRank foreignOverlapByRank_;
    
    // stores the overlap the current processes has in the domain
    // "owned" by some remote ranks
    DomesticOverlap domesticOverlap_;

    // stores the "front" indices of the domestic overlap
    DomesticFront domesticOverlapFront_;

    // stores the "border" indices of the domestic overlap
    DomesticBorder domesticOverlapBorder_;

    // number of domestic local indices
    int numDomesticIndices_;
    
    // offset for the overlap indices for each peer process in number
    // of indices after the last "owned" index
    std::map<ProcessRank, int> domesticOffset_;

    // maps a remote local index to a domestic local index
    std::map<ProcessRank, std::map<RowIndex, RowIndex> > foreignToLocalMap_;
};

} // namespace Dumux

#endif
