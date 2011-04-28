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
    
        assign(M);
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
        sync_();
    }

    void print() const 
    {
        overlap_.print();
        
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
        Dune::printmatrix(std::cout, *static_cast<const BCRSMatrix*>(this), "M", "row");
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

    bool isAdditionalEntry_(bool isPeerLocalRow,
                            int peerRank,
                            int colIdx)
    {
        // indices for which the current process is not master are
        // never additional indices for the peer (only the master
        // processes sends such entries)
        if (!overlap_.iAmMasterOf(colIdx))
            return false;
        
        if (isPeerLocalRow) {
            // for rows which do already exist at the peer rank
            // because they are border rows, an additional index
            // is required only if the column index is an
            // artificial one
            if (!overlap_.isRemoteIndexFor(peerRank, colIdx)) 
                return false;
        }
        else {
            // for rows which do not already exist at the peer
            // rank, all columns which are seen by the peer need
            // to be added.
            if (!overlap_.isDomesticIndexFor(peerRank, colIdx)) 
                return false;
        }
        
        // all other entries are additional
        return true;
    }

    int numAdditionalEntriesInRow_(const BCRSMatrix &M, 
                                   int peerRank, 
                                   int rowIdx)
    {
        bool isPeerLocalRow = overlap_.isLocalIndexFor(peerRank, rowIdx);

        int numEntries = 0;
        typedef typename BCRSMatrix::ConstColIterator ColIt;
        ColIt colIt = M[rowIdx].begin();
        ColIt colEndIt = M[rowIdx].end();
        for (; colIt != colEndIt; ++colIt) {
            if (isAdditionalEntry_(isPeerLocalRow, peerRank, colIt.index()))
                ++numEntries;
        }

        return numEntries;

    };

    int numDomesticEntriesInRowFor_(const BCRSMatrix &M, int peerRank, int rowIdx)
    {
        int numEntries = 0;

        typedef typename BCRSMatrix::ConstColIterator ColIt;
        ColIt colIt = M[rowIdx].begin();
        ColIt colEndIt = M[rowIdx].end();
        for (; colIt != colEndIt; ++colIt) {
            if (overlap_.isDomesticIndexFor(peerRank, colIt.index()))
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
            int numEntries = numAdditionalEntriesInRow_(M, peerRank, rowIdx);

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
            bool isLocalRowForPeer = overlap_.isLocalIndexFor(peerRank, rowIdx);
            
            // loop over the columns of the matrix' row and find out
            // the non-border entries
            int numEntries = numAdditionalEntriesInRow_(M, peerRank, rowIdx);

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
                if (isAdditionalEntry_(isLocalRowForPeer, peerRank, colIt.index()))
                {
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
        }
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
                if (!overlap_.isLocal(domRowIdx) || 
                    !overlap_.isLocal(domColIdx))
                {
                    this->addindex(domRowIdx, domColIdx);
                }
            }

            delete[] recvBuff2;
        };
    
    };

    // communicates and adds up the contents of overlapping rows
    void sync_()
    {
        // first, send all entries to the peers with lower ranks
        typename PeerSet::const_iterator peerIt = overlap_.peerSet().begin();
        typename PeerSet::const_iterator peerEndIt = overlap_.peerSet().end();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            
            if (peerRank < myRank_)
                sendEntries_(peerRank);
        }

        // then recieve entries from the peers with higher ranks
        peerIt = overlap_.peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            
            if (peerRank > myRank_)
                receiveEntries_(peerRank);
        }

        // then send all entries to the peers with higher ranks
        peerIt = overlap_.peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            
            if (peerRank > myRank_)
                sendEntries_(peerRank);
        }

        // finally, receive all entries from peers with lower ranks
        peerIt = overlap_.peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            
            if (peerRank < myRank_)
                receiveEntries_(peerRank);
        }
    }

    void sendEntries_(int peerRank)
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
        
        std::cout << "send " << numOverlapRows << " rows to " << peerRank << "\n";
        typename ForeignOverlapWithPeer::const_iterator it = peerOverlap.begin();
        typename ForeignOverlapWithPeer::const_iterator endIt = peerOverlap.end();
        for (; it != endIt; ++it) {
            int rowIdx = std::get<0>(*it);
            
            // loop over the columns of the matrix' row and find out
            // the number of entries which the peer sees
            int numEntries = numDomesticEntriesInRowFor_(*this, peerRank, rowIdx);

            // send the (global row index, number of entries) pair to
            // the peer
            int sendBuff[2] = { overlap_.domesticToGlobal(rowIdx), numEntries };
            MPI_Send(sendBuff, // buff
                     2, // count
                     MPI_INT, // data type
                     peerRank, 
                     0, // tag
                     MPI_COMM_WORLD); // communicator

            std::cout << "   send " << numEntries << " entries to " << peerRank << " for row " << sendBuff[0] << "\n";

            typedef typename BCRSMatrix::block_type MatrixBlock;
            int *indicesSendBuff = new int[numEntries];
            MatrixBlock *valuesSendBuff = new MatrixBlock[numEntries];
            typedef typename BCRSMatrix::ConstColIterator ColIt;
            int i = 0;
            ColIt colIt = (*this)[rowIdx].begin();
            ColIt colEndIt = (*this)[rowIdx].end();
            for (; colIt != colEndIt; ++colIt) {
                if (overlap_.isDomesticIndexFor(peerRank, colIt.index())) {
                    indicesSendBuff[i] = overlap_.domesticToGlobal(colIt.index());
                    std::cout << "       send column " << indicesSendBuff[i] << "\n";
                    i += 1;
                }
            }

            MPI_Send(indicesSendBuff, // buff
                     numEntries, // count
                     MPI_INT, // data type
                     peerRank, 
                     0, // tag
                     MPI_COMM_WORLD); // communicator

            MPI_Send(valuesSendBuff, // buff
                     numEntries * sizeof(MatrixBlock), // count
                     MPI_BYTE, // data type
                     peerRank,
                     0, // tag
                     MPI_COMM_WORLD); // communicator

            delete[] indicesSendBuff;
            delete[] valuesSendBuff;
        };
    }

    void receiveEntries_(int peerRank)
    {
        // receive size of foreign overlap to peer
        int numOverlapRows;
        MPI_Recv(&numOverlapRows, // buff
                 1, // count
                 MPI_INT, // data type
                 peerRank, 
                 0, // tag
                 MPI_COMM_WORLD, // communicator
                 MPI_STATUS_IGNORE);
        std::cout << "receive " << numOverlapRows << " rows from " << peerRank << "\n";
        
        for (int i = 0; i < numOverlapRows; ++i) {
            // send the (global row index, number of entries) pair to
            // the peer
            int recvBuff[2];
            MPI_Recv(recvBuff, // buff
                     2, // count
                     MPI_INT, // data type
                     peerRank, 
                     0, // tag
                     MPI_COMM_WORLD, // communicator
                     MPI_STATUS_IGNORE);

            // loop over the columns of the matrix' row and find out
            // the number of entries which the peer sees
            int domRowIdx = overlap_.globalToDomestic(recvBuff[0]);
            int numEntries = recvBuff[1];
            std::cout << "   receive " << numEntries << " entries from " << peerRank << " for row " << recvBuff[0] << "\n";

            typedef typename BCRSMatrix::block_type MatrixBlock;
            int *indicesRecvBuff = new int[numEntries];
            MatrixBlock *valuesRecvBuff = new MatrixBlock[numEntries];
            MPI_Recv(indicesRecvBuff, // buff
                     numEntries, // count
                     MPI_INT, // data type
                     peerRank, 
                     0, // tag
                     MPI_COMM_WORLD, // communicator
                     MPI_STATUS_IGNORE);

            MPI_Recv(valuesRecvBuff, // buff
                     numEntries * sizeof(MatrixBlock), // count
                     MPI_BYTE, // data type
                     peerRank,
                     0, // tag
                     MPI_COMM_WORLD,// communicator
                     MPI_STATUS_IGNORE); 
            
            for (int j = 0; j < numEntries; ++j) {
                int domColIdx = overlap_.globalToDomestic(indicesRecvBuff[j]);
                std::cout << "       receive column " << indicesRecvBuff[j] << "\n";
                std::cout.flush();

                if (overlap_.iAmMasterOf(domRowIdx))
                    (*this)[domRowIdx][domColIdx] += valuesRecvBuff[j];
                else
                    (*this)[domRowIdx][domColIdx] = valuesRecvBuff[j];
            }

            delete[] indicesRecvBuff;
            delete[] valuesRecvBuff;
        };
    }

    int myRank_;
    Overlap overlap_;
    std::map<ProcessRank, int> numOverlapEntries_;
    OverlapStructure overlapEntries_;
};

} // namespace Dumux

#endif
