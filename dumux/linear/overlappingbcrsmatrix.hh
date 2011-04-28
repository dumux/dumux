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
#ifndef DUMUX_OVERLAPPING_BCRS_MATRIX_HH
#define DUMUX_OVERLAPPING_BCRS_MATRIX_HH

#include <dune/istl/scalarproducts.hh>

#include <algorithm>
#include <list>
#include <set>
#include <map>
#include <memory>

namespace Dumux {

template <class BCRSMatrix>
class OverlappingBCRSMatrix : public BCRSMatrix
{
    typedef BCRSMatrix ParentType;

public:
    typedef Dumux::DomesticOverlapFromBCRSMatrix<BCRSMatrix> Overlap;
    
private:
    typedef typename Overlap::Index RowIndex;
    typedef typename Overlap::Index ColIndex;
    typedef typename Overlap::PeerSet PeerSet;
    typedef typename Overlap::BorderList BorderList;
    typedef typename Overlap::ProcessRank ProcessRank;
    typedef typename Overlap::ForeignOverlapWithPeer ForeignOverlapWithPeer;
    typedef std::vector<ColIndex> PeerColumns;
    typedef std::pair<RowIndex,  PeerColumns> PeerRow;
    typedef std::vector<PeerRow> PeerRows;

    typedef std::vector<std::set<ColIndex> > Entries;

public:
    typedef typename ParentType::RowIterator RowIterator;
    typedef typename ParentType::ColIterator ColIterator;
    typedef typename ParentType::ConstColIterator ConstColIterator;
    typedef typename ParentType::field_type field_type;
    typedef typename ParentType::block_type block_type;

    OverlappingBCRSMatrix(const BCRSMatrix &M, 
                          const BorderList &borderList, 
                          int overlapSize)
    {
        overlap_ = std::shared_ptr<Overlap>(new Overlap(M, borderList, overlapSize));
        overlap_->print();

        MPI_Comm_rank(MPI_COMM_WORLD, &myRank_);

        // build the overlapping matrix from the non-overlapping
        // matrix and the overlap
        build_(M);
    
        assign(M);
    };
    
    OverlappingBCRSMatrix(const OverlappingBCRSMatrix &M)
        : ParentType(M)
        , overlap_(M.overlap_)
    {
    };

    /*!
     * \brief Returns the domestic overlap for the process.
     */
    const Overlap &overlap() const
    { return *overlap_; }

    /*!
     * \brief Assign and syncronize the overlapping matrix from a
     *       non-overlapping one.
     */
    void assign(const BCRSMatrix &M)
    {
        // first, set everything to 0
        BCRSMatrix::operator=(0.0);

        // assign the local rows
        for (int rowIdx = 0; rowIdx < M.N(); ++rowIdx) {
            ConstColIterator colIt = M[rowIdx].begin();
            ConstColIterator colEndIt = M[rowIdx].end();
            ColIterator myColIt = (*this)[rowIdx].begin();
            for (; colIt != colEndIt; ++colIt) {
                while (true) {
                    if (myColIt.index() == colIt.index())
                        break;

                    ++ myColIt;
                }
                assert(myColIt.index() == colIt.index());
                
                (*myColIt) = *colIt;
            }
        }
        
        // communicate and add the contents of overlapping rows
        sync_();
    }

    void print() const 
    {
        overlap_->print();
        
        for (int i = 0; i < this->N(); ++i) {
            if (overlap_->isLocal(i))
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
                    if (overlap_->isBorder(j))
                        std::cout << "|";
                    else if (overlap_->isLocal(j))
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
        int numDomestic = overlap_->numDomestic();
                       
        // allocate the rows
        this->setSize(numDomestic, numDomestic);
        this->setBuildMode(ParentType::random);
       
        // communicate the entries
        buildIndices_(M);
    }

    int numDomesticEntriesInRowFor_(const BCRSMatrix &M, int peerRank, int rowIdx)
    {
        int numEntries = 0;

        typedef typename BCRSMatrix::ConstColIterator ColIt;
        ColIt colIt = M[rowIdx].begin();
        ColIt colEndIt = M[rowIdx].end();
        for (; colIt != colEndIt; ++colIt) {
            if (overlap_->isDomesticIndexFor(peerRank, colIt.index()))
                ++numEntries;
        }

        return numEntries;
    }

    void buildIndices_(const BCRSMatrix &M)
    {
        /////////
        // first, add all local matrix entries
        /////////
        entries_.resize(overlap_->numDomestic());
        for (int rowIdx = 0; rowIdx < M.N(); ++rowIdx) {
            ConstColIterator colIt = M[rowIdx].begin();
            ConstColIterator colEndIt = M[rowIdx].end();
            for (; colIt != colEndIt; ++colIt) {
                entries_[rowIdx].insert(colIt.index());
            }
        }
        
        /////////
        // add the indices for all additional entries
        /////////

        // first, send all indices to the peers with lower ranks
        typename PeerSet::const_iterator peerIt = overlap_->peerSet().begin();
        typename PeerSet::const_iterator peerEndIt = overlap_->peerSet().end();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            
            if (peerRank < myRank_)
                sendRowIndices_(M, peerRank);
        }

        // then recieve the indices from the peers with higher ranks
        peerIt = overlap_->peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            
            if (peerRank > myRank_)
                receiveRowIndices_(peerRank);
        }

        // then send all indices to the peers with higher ranks
        peerIt = overlap_->peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            
            if (peerRank > myRank_)
                sendRowIndices_(M, peerRank);
        }

        // then receive all indices from peers with lower ranks
        peerIt = overlap_->peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            
            if (peerRank < myRank_)
                receiveRowIndices_(peerRank);
        }
        
        /////////
        // actually initialize the BCRS matrix structure
        /////////

        // set the row sizes
        int numDomestic = overlap_->numDomestic();
        for (int rowIdx = 0; rowIdx < numDomestic; ++rowIdx) {
            this->setrowsize(rowIdx, entries_[rowIdx].size());
        }
        this->endrowsizes();

        // set the indices
        for (int rowIdx = 0; rowIdx < numDomestic; ++rowIdx) {
            const std::set<ColIndex> &colIndices = entries_[rowIdx];

            typename std::set<ColIndex>::const_iterator colIdxIt = colIndices.begin();
            typename std::set<ColIndex>::const_iterator colIdxEndIt = colIndices.end();
            for (; colIdxIt != colIdxEndIt; ++colIdxIt)
                this->addindex(rowIdx, *colIdxIt);
        }
        this->endindices();

        // free the memory occupied by the array of the matrix entries
        entries_.resize(0);
    }

    // send the overlap indices to a peer
    void sendRowIndices_(const BCRSMatrix &M, int peerRank)
    {
        // send the number of non-border entries in the matrix
        const ForeignOverlapWithPeer &peerOverlap = overlap_->foreignOverlapWithPeer(peerRank);

        // send size of foreign overlap to peer
        int numOverlapRows = peerOverlap.size();
        MPI_Bsend(&numOverlapRows, // buff
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
            int numEntries = numDomesticEntriesInRowFor_(M, peerRank, rowIdx);

            // send the (global row index, number of additional
            // entries) pair to the peer
            int sendBuff[2] = { overlap_->domesticToGlobal(rowIdx), numEntries };
            MPI_Bsend(sendBuff, // buff
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
                if (overlap_->isDomesticIndexFor(peerRank, rowIdx) && 
                    overlap_->isDomesticIndexFor(peerRank, colIt.index()))
                {
                    sendBuff2[i] = overlap_->domesticToGlobal(colIt.index());
                    i += 1;
                }
            }
        
            MPI_Bsend(sendBuff2, // buff
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

            int domRowIdx = overlap_->globalToDomestic(recvBuff[0]);
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
                int domColIdx = overlap_->globalToDomestic(recvBuff2[i]);
                entries_[domRowIdx].insert(domColIdx);
            }

            delete[] recvBuff2;
        };
    
    };

    // communicates and adds up the contents of overlapping rows
    void sync_()
    {
        // first, recieve entries from the peers with higher ranks
        typename PeerSet::const_iterator peerIt = overlap_->peerSet().begin();
        typename PeerSet::const_iterator peerEndIt = overlap_->peerSet().end();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            
            if (peerRank > myRank_)
                receiveEntries_(peerRank);
        }

        // then, send all entries to the peers with lower ranks
        peerIt = overlap_->peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            
            if (peerRank < myRank_)
                sendEntries_(peerRank);
        }

        // then, receive all entries from peers with lower ranks
        peerIt = overlap_->peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            
            if (peerRank < myRank_)
                receiveEntries_(peerRank);
        }

        // finally, send all entries to the peers with higher ranks
        peerIt = overlap_->peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            
            if (peerRank > myRank_)
                sendEntries_(peerRank);
        }
    }

    void sendEntries_(int peerRank)
    {
        // send the number of non-border entries in the matrix
        const ForeignOverlapWithPeer &peerOverlap = overlap_->foreignOverlapWithPeer(peerRank);

        // send size of foreign overlap to peer
        int numOverlapRows = peerOverlap.size();
        MPI_Bsend(&numOverlapRows, // buff
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
            // the number of entries which the peer sees
            int numEntries = numDomesticEntriesInRowFor_(*this, peerRank, rowIdx);

            // send the (global row index, number of entries) pair to
            // the peer
            int sendBuff[2] = { overlap_->domesticToGlobal(rowIdx), numEntries };
            MPI_Bsend(sendBuff, // buff
                     2, // count
                     MPI_INT, // data type
                     peerRank, 
                     0, // tag
                     MPI_COMM_WORLD); // communicator

            typedef typename BCRSMatrix::block_type MatrixBlock;
            int *indicesSendBuff = new int[numEntries];
            MatrixBlock *valuesSendBuff = new MatrixBlock[numEntries];
            typedef typename BCRSMatrix::ConstColIterator ColIt;

            int i = 0;
            ColIt colIt = (*this)[rowIdx].begin();
            ColIt colEndIt = (*this)[rowIdx].end();
            for (; colIt != colEndIt; ++colIt) {
                if (overlap_->isDomesticIndexFor(peerRank, colIt.index())) {
                    indicesSendBuff[i] = overlap_->domesticToGlobal(colIt.index());
                    valuesSendBuff[i] = (*colIt);
                    i += 1;
                }
            }

            MPI_Bsend(indicesSendBuff, // buff
                     numEntries, // count
                     MPI_INT, // data type
                     peerRank, 
                     0, // tag
                     MPI_COMM_WORLD); // communicator

            MPI_Bsend(valuesSendBuff, // buff
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
            int domRowIdx = overlap_->globalToDomestic(recvBuff[0]);
            int numEntries = recvBuff[1];

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
                int domColIdx = overlap_->globalToDomestic(indicesRecvBuff[j]);

                if (peerRank < myRank_)
                    (*this)[domRowIdx][domColIdx] = valuesRecvBuff[j];
                else
                    (*this)[domRowIdx][domColIdx] += valuesRecvBuff[j];
            }
            

            delete[] indicesRecvBuff;
            delete[] valuesRecvBuff;
        };
    }

    int myRank_;
    std::shared_ptr<Overlap> overlap_;
    std::map<ProcessRank, int> numOverlapEntries_;
    Entries entries_;
    
};

} // namespace Dumux

#endif
