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
#ifndef DUMUX_OVERLAPPING_MATRIX_HH
#define DUMUX_OVERLAPPING_MATRIX_HH

#include <dune/istl/scalarproducts.hh>

#include <algorithm>
#include <list>
#include <set>
#include <map>

namespace Dumux {

template <class BCRSMatrix>
class OverlappingBCRSMatrix : public BCRSMatrix
{
    typedef BCRSMatrix ParentType;
    typedef Dumux::DomesticOverlapFromBCRSMatrix<BCRSMatrix> Overlap;

    typedef typename Overlap::Index RowIndex;
    typedef typename Overlap::Index ColIndex;
    typedef typename Overlap::PeerSet PeerSet;
    typedef typename Overlap::BorderList BorderList;
    typedef typename Overlap::ProcessRank ProcessRank;
    typedef typename Overlap::ForeignOverlapWithPeer ForeignOverlapWithPeer;
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

    OverlappingBCRSMatrix(const BCRSMatrix &M, 
                          const BorderList &borderList, 
                          int overlapSize)
        : overlap_(M, borderList, overlapSize)
    {
        overlap_.print();

        MPI_Comm_rank(MPI_COMM_WORLD, &myRank_);

        // build the overlapping matrix from the non-overlapping
        // matrix and the overlap
        build_(M);
    
        //sync(M);
    };
    
    OverlappingBCRSMatrix(const OverlappingBCRSMatrix &M)
        : ParentType(M), 
          overlap_(M.overlap_)
    {
    };

    /*!
     * \brief Assign and syncronize the overlapping matrix from a
     *       non-overlapping one.
     */
    void assign(const BCRSMatrix &M)
    {
        // set everything to 0
        ParentType::operator=(0);

        // assign the local rows
        for (int rowIdx = 0; rowIdx < M.N(); ++rowIdx) {
            ConstColIterator colIt = M[rowIdx].begin();
            ConstColIterator colEndIt = M[rowIdx].end();
            ColIterator myColIt = (*this)[rowIdx].begin();
            for (; colIt != colEndIt; ++colIt) {
                while (myColIt.index() < colIt.index()) {
                    ++ myColIt;
                }
                assert(myColIt.index() == colIt.index());
                
                (*myColIt) = *colIt;
            }
        };
        
        // communicate and add the contents of overlapping rows
        sync_(M);
    }

    void print() const 
    {
        for (int i = 0; i < this->N(); ++i) {
            if (overlap_.isLocal(i))
                std::cout << " ";
            else 
                std::cout << "*";
            std::cout << "row " << i << " ";
            
            typedef typename BCRSMatrix::ConstColIterator ColIt;
            ColIt colIt = (*this)[i].begin();
            ColIt colEndIt = (*this)[i].end();
            for (int j = 0; j < this->M(); ++j) {
                if (colIt != colEndIt && j == colIt.index()) { 
                    ++colIt;
                    if (overlap_.isBorder(j))
                        std::cout << "|";
                    else if (overlap_.isLocal(j))
                        std::cout << "X";
                    else
                        std::cout << "*";
                }
                else
                    std::cout << " ";
            }
            std::cout << "\n";
        };
    };

private:
    void build_(const BCRSMatrix &M)
    {
        int numDomestic = overlap_.numDomestic();
                       
        // allocate the rows
        this->setSize(numDomestic, numDomestic);
        this->setBuildMode(ParentType::random);
       
        // set the row size for all additional rows
        buildRowSizes_(M);

        // communicate the entries
        buildIndices_(M);
    }

    void buildRowSizes_(const BCRSMatrix &M)
    {
        // copy the rows for the local indices
        int numLocal = M.N();
        for (int rowIdx = 0; rowIdx < numLocal; ++rowIdx)
            this->setrowsize(rowIdx, M.getrowsize(rowIdx));

        // initialize the remaining row sizes to 0
        int numDomestic = overlap_.numDomestic();
        for (int rowIdx = numLocal; rowIdx < numDomestic; ++rowIdx)
            this->setrowsize(rowIdx, 0);

        typename PeerSet::const_iterator peerIt = overlap_.peerSet().begin();
        typename PeerSet::const_iterator peerEndIt = overlap_.peerSet().end();

        // send the number of additional entries for each row to the
        // peers with lower ranks
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            
            if (peerRank < myRank_)
                sendRowSizesToPeer_(M, peerRank);
        }

        // receive the number of additional entries for each row from
        // the peers with higher ranks
        peerIt = overlap_.peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            
            if (peerRank > myRank_)
                receiveRowSizesFromPeer_(peerRank);
        }

        // send the number of additional entries for each row to
        // the peers with higher ranks
        peerIt = overlap_.peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            
            if (peerRank > myRank_)
                sendRowSizesToPeer_(M, peerRank);
        }

        // receive the number of additional entries for each row from
        // the peers with lower ranks
        peerIt = overlap_.peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            
            if (peerRank < myRank_)
                receiveRowSizesFromPeer_(peerRank);
        }
        this->endrowsizes();
    }

    int numRemoteEntriesInRow_(const BCRSMatrix &M, int peerRank, int rowIdx)
    {
        int numEntries = 0;

        typedef typename BCRSMatrix::ConstColIterator ColIt;
        ColIt colIt = M[rowIdx].begin();
        ColIt colEndIt = M[rowIdx].end();
        for (; colIt != colEndIt; ++colIt) {
            if (overlap_.isRemoteIndexFor(peerRank, colIt.index()))
                ++numEntries;
        }

        return numEntries;
    }

    void sendRowSizesToPeer_(const BCRSMatrix &M, int peerRank)
    {
        // send the number of non-border entries in the matrix
        const ForeignOverlapWithPeer &peerOverlap = overlap_.foreignOverlapWithPeer(peerRank);

        // send size of foreign overlap to peer
        int numOverlapRows = peerOverlap.size();
        MPI_Send(&numOverlapRows, // buff
                 1, // count
                 MPI_INT, // data type
                 peerRank, 
                 0, // tag
                 MPI_COMM_WORLD); // communicator
        

        typename ForeignOverlapWithPeer::const_iterator it = peerOverlap.begin();
        typename ForeignOverlapWithPeer::const_iterator endIt = peerOverlap.end();
        for (; it != endIt; ++it) {
            int rowIdx = std::get<0>(*it);
            
            // loop over the columns of the matrix' row and find out
            // the non-border entries
            int numEntries = numRemoteEntriesInRow_(M, peerRank, rowIdx);

            // send the (global row index, number of additional
            // entries) pair to the peer
            int sendBuff[2] = { overlap_.domesticToGlobal(rowIdx), numEntries };
            MPI_Send(sendBuff, // buff
                     2, // count
                     MPI_INT, // data type
                     peerRank, 
                     0, // tag
                     MPI_COMM_WORLD); // communicator
        };
    }

    void receiveRowSizesFromPeer_(int peerRank)
    {
        // send size of foreign overlap to peer
        int numOverlapRows;
        MPI_Recv(&numOverlapRows, // buff
                 1, // count
                 MPI_INT, // data type
                 peerRank, 
                 0, // tag
                 MPI_COMM_WORLD, // communicator
                 MPI_STATUS_IGNORE); 
       
        for (int i = 0; i < numOverlapRows; ++i) {
            // receive the (global row index, number of additional
            // entries) pair from the peer
            int recvBuff[2];
            MPI_Recv(recvBuff, // buff
                     2, // count
                     MPI_INT, // data type
                     peerRank, 
                     0, // tag
                     MPI_COMM_WORLD, // communicator
                     MPI_STATUS_IGNORE);

            // translate global index to local one and increment row
            // size
            int globalIdx = recvBuff[0];
            int domesticIdx = overlap_.globalToDomestic(globalIdx);
            int numEntries = this->getrowsize(domesticIdx) + recvBuff[1];

            this->setrowsize(domesticIdx, numEntries);
        };
    }

    void buildIndices_(const BCRSMatrix &M)
    {
        // add the indices for the local entries
        // copy the rows for the local indices
        int numLocal = M.N();
        for (int rowIdx = 0; rowIdx < numLocal; ++rowIdx) {
            ConstColIterator colIt = M[rowIdx].begin();
            ConstColIterator colEndIt = M[rowIdx].end();
            for (; colIt != colEndIt; ++colIt)
                this->addindex(rowIdx, colIt.index());
        }

        // add the indices for all additional entries

        // first, send all indices to the peers with lower ranks
        typename PeerSet::const_iterator peerIt = overlap_.peerSet().begin();
        typename PeerSet::const_iterator peerEndIt = overlap_.peerSet().end();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            
            if (peerRank < myRank_)
                sendRowIndices_(M, peerRank);
        }

        // then recieve the indices from the peers with higher ranks
        peerIt = overlap_.peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            
            if (peerRank > myRank_)
                receiveRowIndices_(peerRank);
        }

        // then send all indices to the peers with higher ranks
        peerIt = overlap_.peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            
            if (peerRank > myRank_)
                sendRowIndices_(M, peerRank);
        }

        // then receive all indices from peers with lower ranks
        peerIt = overlap_.peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            
            if (peerRank < myRank_)
                receiveRowIndices_(peerRank);
        }

        this->endindices();
    }

    // send the overlap indices to a peer
    void sendRowIndices_(const BCRSMatrix &M, int peerRank)
    {
        // send the number of non-border entries in the matrix
        const ForeignOverlapWithPeer &peerOverlap = overlap_.foreignOverlapWithPeer(peerRank);

        // send size of foreign overlap to peer
        int numOverlapRows = peerOverlap.size();
        MPI_Send(&numOverlapRows, // buff
                 1, // count
                 MPI_INT, // data type
                 peerRank, 
                 0, // tag
                 MPI_COMM_WORLD); // communicator
        
        typename ForeignOverlapWithPeer::const_iterator it = peerOverlap.begin();
        typename ForeignOverlapWithPeer::const_iterator endIt = peerOverlap.end();
        for (; it != endIt; ++it) {
            int rowIdx = std::get<0>(*it);
            
            // loop over the columns of the matrix' row and find out
            // the non-border entries
            int numEntries = numRemoteEntriesInRow_(M, peerRank, rowIdx);

            // send the (global row index, number of additional
            // entries) pair to the peer
            int sendBuff[2] = { overlap_.domesticToGlobal(rowIdx), numEntries };
            MPI_Send(sendBuff, // buff
                     2, // count
                     MPI_INT, // data type
                     peerRank, 
                     0, // tag
                     MPI_COMM_WORLD); // communicator

            int *sendBuff2 = new int[numEntries];
            typedef typename BCRSMatrix::ConstColIterator ColIt;
            int i = 0;
            ColIt colIt = M[rowIdx].begin();
            ColIt colEndIt = M[rowIdx].end();
            for (; colIt != colEndIt; ++colIt) {
                if (overlap_.isRemoteIndexFor(peerRank, colIt.index())) {
                    sendBuff2[i] = overlap_.domesticToGlobal(colIt.index());
                    i += 1;
                }
            }

            MPI_Send(sendBuff2, // buff
                     numEntries, // count
                     MPI_INT, // data type
                     peerRank, 
                     0, // tag
                     MPI_COMM_WORLD); // communicator

            delete[] sendBuff2;
        };
    
    };

    // receive the overlap indices to a peer
    void receiveRowIndices_(int peerRank)
    {
        // receive the of the domestic overlap to peer
        int numOverlapRows;
        MPI_Recv(&numOverlapRows, // buff
                 1, // count
                 MPI_INT, // data type
                 peerRank, 
                 0, // tag
                 MPI_COMM_WORLD, // communicator
                 MPI_STATUS_IGNORE);
        for (int j = 0; j < numOverlapRows; ++j) {
            // receive the (global row index, number of additional
            // entries) pair from the peer
            int recvBuff[2];
            MPI_Recv(recvBuff, // buff
                     2, // count
                     MPI_INT, // data type
                     peerRank, 
                     0, // tag
                     MPI_COMM_WORLD, // communicator
                     MPI_STATUS_IGNORE);

            int domRowIdx = overlap_.globalToDomestic(recvBuff[0]);
            int numEntries = recvBuff[1];
            
            // receive the actual column indices
            int *recvBuff2 = new int[numEntries];
            MPI_Recv(recvBuff2, // buff
                     numEntries, // count
                     MPI_INT, // data type
                     peerRank, 
                     0, // tag
                     MPI_COMM_WORLD, // communicator
                     MPI_STATUS_IGNORE);
            for (int i = 0; i < numEntries; ++i) {
                int domColIdx = overlap_.globalToDomestic(recvBuff2[i]);
                this->addindex(domRowIdx, domColIdx);
            }

            delete[] recvBuff2;
        };
    
    };

    // communicates and adds up the contents of overlapping rows
    void sync_()
    {
    }

/*
    // retrieve the number of additional entries for each row in the
    // domestic overlap. this communicates via MPI
    void communicateOverlapStructure_(const BCRSMatrix &M)
    {
        const PeerSet &peerSet = overlap_.peerSet();
        
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
            overlap_.foreignOverlapWithPeer(peerRank).size();
        int *buff = new int[msgSize];
        const ForeignOverlapWithPeer &olist = overlap_.foreignOverlapWithPeer(peerRank);
        typename ForeignOverlapWithPeer::const_iterator it = olist.begin();
        typename ForeignOverlapWithPeer::const_iterator endIt = olist.end();
        for (int i = 0; it != endIt; ++it, ++i) {
            int rowIdx = std::get<0>(*it);
            buff[2*i + 0] = rowIdx;
            if (overlap_.isForeignFront(peerRank, rowIdx))
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

        const ForeignOverlapWithPeer &olist = overlap_.foreignOverlapWithPeer(peerRank);
        
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
            
            if (overlap_.isForeignFront(peerRank, rowIdx)) {
                continue; // we do not send anything for front rows
            }

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
        int numRows = overlap_.domesticOverlapWithPeer(peerRank).size();
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

            if (overlap_.isFront(peerRank, rowIdx))
                overlapEntries_[peerRank][i].second.resize(0);
        }
        delete[] buff;
    };
    
    void syncronizeEntryValues_(const BCRSMatrix &M)
    {
        const PeerSet &peerSet = overlap_.peerSet();
        
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

            MPI_Isend(sendMsgBuff[i], // pointer to user data 
                      msgSize*sizeof(field_type), // size of user data array
                      MPI_BYTE, // type of user data
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
*/

/*
    field_type *createValuesMsg_(const BCRSMatrix &M, int &msgSize, int peerRank)
    {
        msgSize = 0;

        const ForeignOverlapWithPeer &olist = overlap_.foreignOverlapWithPeer(peerRank);
        
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
            
            if (overlap_.isForeignFront(peerRank, rowIdx))
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
        const DomesticOverlapWithPeer &olist = overlap_.domesticOverlapWithPeer(peerRank);
        int numRows = olist.size();
        int msgSize = numOverlapEntries_[peerRank] * block_type::rows*block_type::cols;
        field_type *buff = new field_type[msgSize];

        MPI_Recv(buff, // receive message buffer
                 msgSize*sizeof(field_type), // receive message size
                 MPI_BYTE, // object type
                 peerRank, // peer rank
                 0, // identifier
                 MPI_COMM_WORLD, // communicator
                 MPI_STATUS_IGNORE); // status

        int pos = 0;
        for (int i = 0; i < numRows; ++i) {
            int rowIdx = overlapEntries_[peerRank][i].first;
            int domRowIdx = overlap_.foreignToDomesticIndex(peerRank, rowIdx);

            if (overlap_.isFront(peerRank, rowIdx)) {
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
                int domColIdx = overlap_.foreignToDomesticIndex(peerRank, colIdx);
                     
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
*/

    int myRank_;
    Overlap overlap_;
    std::map<ProcessRank, int> numOverlapEntries_;
    OverlapStructure overlapEntries_;
};

} // namespace Dumux

#endif
