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
 * \brief A block vector which creates an algebraic overlap of
 *        arbitrary size.
 */
#ifndef DUMUX_OVERLAPPING_BLOCK_VECTOR_HH
#define DUMUX_OVERLAPPING_BLOCK_VECTOR_HH

#include <dune/istl/bvector.hh>

namespace Dumux {

template <class FieldVector, class Overlap>
class OverlappingBlockVector
    : public Dune::BlockVector<FieldVector>
{
    typedef Dune::BlockVector<FieldVector> ParentType;
    typedef Dune::BlockVector<FieldVector> BlockVector;

    typedef typename Overlap::Index RowIndex;
    typedef typename Overlap::ProcessRank ProcessRank;
    typedef typename Overlap::PeerSet PeerSet;
    typedef typename Overlap::ForeignOverlapWithPeer ForeignOverlapWithPeer;
    typedef typename Overlap::DomesticOverlapWithPeer DomesticOverlapWithPeer;

    typedef typename ParentType::field_type Scalar;

public:
    OverlappingBlockVector(const BlockVector &nbv, const Overlap &overlap)
        : ParentType(overlap.numDomestic())
        , overlap_(&overlap)
    {
        createSendBuffers_();
        assignAdd(nbv);
    };

    /*!
     * \brief Copy constructor.
     */
    OverlappingBlockVector(const OverlappingBlockVector &obv)
        : ParentType(obv)
        , overlap_(obv.overlap_)
    {
        createSendBuffers_();
    }

    ~OverlappingBlockVector()
    {
        // delete buffer pointers
        typename PeerSet::const_iterator peerIt;
        typename PeerSet::const_iterator peerEndIt = overlap_->peerSet().end();
        
        // send all entries to all peers
        peerIt = overlap_->peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;

            delete[] indicesSendBuff_[peerRank];
            delete[] valuesSendBuff_[peerRank];
        }
    }
    
    /*!
     * \brief Assign an overlapping block vector from a
     *        non-overlapping one, border entries are added.
     */
    void assignAdd(const BlockVector &nbv)
    {
        // assign the local rows
        int numLocal = overlap_->numLocal();
        for (int rowIdx = 0; rowIdx < numLocal; ++rowIdx) {
            (*this)[rowIdx] = nbv[rowIdx];
        }

        int numDomestic = overlap_->numDomestic();
        for (int rowIdx = numLocal; rowIdx < numDomestic; ++rowIdx) {
            (*this)[rowIdx] = 0;
        }
        
        // add up the contents of overlapping rows
        syncAdd();
    }

    /*!
     * \brief Assign the local values to a non-overlapping block
     *        vector.
     */
    void assignTo(BlockVector &nbv) const
    {
        // assign the local rows
        int numLocal = overlap_->numLocal();
        for (int rowIdx = 0; rowIdx < numLocal; ++rowIdx) {
            nbv[rowIdx] = (*this)[rowIdx];
        }
    }

    /*!
     * \brief Syncronize an overlapping block vector by adding up all
     *        overlapping entries.
     */
    void syncAdd()
    {
        typename PeerSet::const_iterator peerIt;
        typename PeerSet::const_iterator peerEndIt = overlap_->peerSet().end();

        // send all entries to all peers
        peerIt = overlap_->peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            sendEntries_(peerRank);
        }

        // recieve all entries to the peers
        peerIt = overlap_->peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            receiveAddEntries_(peerRank);
        }
    };

    using ParentType::operator=;
    /*!
     * \brief Copy constructor.
     */
    OverlappingBlockVector &operator=(const OverlappingBlockVector &obv)
    {
        ParentType::operator=(obv);
        overlap_ = obv.overlap_;
        return *this;
    }
    
    /*!
     * \brief Syncronize an overlapping block vector and take the
     *        arthmetic mean of the entry values of all processes.
     */
    void syncAverage()
    {
        syncAdd();

        int numDomestic = overlap_->numDomestic();
        for (int i = 0; i < numDomestic; ++i) {
            (*this)[i] /= overlap_->numPeers(i) + 1;
        }
    };

    /*!
     * \brief Set all front entities to a given scalar value
     */
    void resetFront(Scalar value = 0.0)
    {
        int numDomestic = this->size();
        for (int i = 0; i < numDomestic; ++i) {
            if (overlap_->isFront(i)) {
                (*this)[i] = value;
            }
        };
    }

    void print() const 
    {
        for (int i = 0; i < this->size(); ++i) {
            std::cout << "row " << i << (overlap_->isLocal(i)?" ":"*") << ": " << (*this)[i] << "\n";
        };
    };

private:
    void createSendBuffers_()
    {
        // create buffer pointers
        typename PeerSet::const_iterator peerIt;
        typename PeerSet::const_iterator peerEndIt = overlap_->peerSet().end();
        
        // send all entries to all peers
        peerIt = overlap_->peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;

            int numEntries = overlap_->domesticOverlapWithPeer(peerRank).size();
            indicesSendBuff_[peerRank] = new RowIndex[numEntries];
            valuesSendBuff_[peerRank] = new FieldVector[numEntries];

            indicesSendReq_[peerRank];
            valuesSendReq_[peerRank];
        }
    }

    void sendEntries_(int peerRank)
    {
        // send the number of non-border entries in the matrix
        const DomesticOverlapWithPeer &domesticOverlap = overlap_->domesticOverlapWithPeer(peerRank);

        // send size of foreign overlap to peer
        int numOverlapRows = domesticOverlap.size();
        MPI_Bsend(&numOverlapRows, // buff
                  1, // count
                  MPI_INT, // data type
                  peerRank, 
                  0, // tag
                  MPI_COMM_WORLD); // communicator
        
        int *indicesSendBuff = indicesSendBuff_[peerRank];
        FieldVector *valuesSendBuff = valuesSendBuff_[peerRank];

        int i = 0;
        typename DomesticOverlapWithPeer::const_iterator domIt = domesticOverlap.begin();
        typename DomesticOverlapWithPeer::const_iterator domEndIt = domesticOverlap.end();
        for (; domIt != domEndIt; ++domIt, ++i) {
            int rowIdx = *domIt;
            indicesSendBuff[i] = overlap_->domesticToGlobal(rowIdx);
            valuesSendBuff[i] = (*this)[rowIdx];
        }

        MPI_Isend(indicesSendBuff, // buff
                  numOverlapRows, // count
                  MPI_INT, // data type
                  peerRank, 
                  0, // tag
                  MPI_COMM_WORLD, // communicator
                  &indicesSendReq_[peerRank]);  // request
        
        MPI_Isend(valuesSendBuff, // buff
                  numOverlapRows * sizeof(FieldVector), // count
                  MPI_BYTE, // data type
                  peerRank,
                  0, // tag
                  MPI_COMM_WORLD,// communicator
                  &valuesSendReq_[peerRank]); // request
    }

    void waitSendFinished_()
    {
        typename PeerSet::const_iterator peerIt;
        typename PeerSet::const_iterator peerEndIt = overlap_->peerSet().end();

        // send all entries to all peers
        peerIt = overlap_->peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            MPI_Wait(&indicesSendReq_[peerRank], MPI_STATUS_IGNORE);
            MPI_Wait(&valuesSendReq_[peerRank], MPI_STATUS_IGNORE);
        }
    }

    void receiveAddEntries_(int peerRank)
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

        int *indicesRecvBuff = new int[numOverlapRows];
        FieldVector *valuesRecvBuff = new FieldVector[numOverlapRows];
        MPI_Recv(indicesRecvBuff, // buff
                 numOverlapRows, // count
                 MPI_INT, // data type
                 peerRank, 
                 0, // tag
                 MPI_COMM_WORLD, // communicator
                 MPI_STATUS_IGNORE);

        MPI_Recv(valuesRecvBuff, // buff
                 numOverlapRows * sizeof(FieldVector), // count
                 MPI_BYTE, // data type
                 peerRank,
                 0, // tag
                 MPI_COMM_WORLD,// communicator
                 MPI_STATUS_IGNORE); 
        
        for (int j = 0; j < numOverlapRows; ++j) {
            int domRowIdx = overlap_->globalToDomestic(indicesRecvBuff[j]);
            (*this)[domRowIdx] += valuesRecvBuff[j];
        }

        delete[] indicesRecvBuff;
        delete[] valuesRecvBuff;
    }

    std::map<ProcessRank, RowIndex *> indicesSendBuff_;
    std::map<ProcessRank, FieldVector *> valuesSendBuff_;

    std::map<ProcessRank, MPI_Request> indicesSendReq_;
    std::map<ProcessRank, MPI_Request> valuesSendReq_;

    const Overlap *overlap_;
};

} // namespace Dumux

#endif
