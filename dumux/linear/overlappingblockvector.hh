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

#include <dumux/parallel/mpibuffer.hh>

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
    OverlappingBlockVector(const Overlap &overlap)
        : ParentType(overlap.numDomestic())
        , overlap_(&overlap)
    {
        createBuffers_();
    };

    /*!
     * \brief Copy constructor.
     */
    OverlappingBlockVector(const OverlappingBlockVector &obv)
        : ParentType(obv)
        , indicesSendBuff_(obv.indicesSendBuff_)
        , indicesRecvBuff_(obv.indicesRecvBuff_)
        , valuesSendBuff_(obv.valuesSendBuff_)
        , valuesRecvBuff_(obv.valuesRecvBuff_)
        , overlap_(obv.overlap_)
    {
    }

    ~OverlappingBlockVector()
    {
    }
    
    /*!
     * \brief Assign an overlapping block vector from a
     *        non-overlapping one, border entries are added.
     */
    void assignAdd(const BlockVector &nbv)
    {
        Valgrind::SetDefined(nbv[7]);

        // assign the local rows
        int numLocal = overlap_->numLocal();
        for (int rowIdx = 0; rowIdx < numLocal; ++rowIdx) {
            (*this)[rowIdx] = nbv[rowIdx];
        }

        // set the remote indices to 0
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
        
        // wait until we have send everything
        waitSendFinished_();
    };

    /*!
     * \brief Syncronize an overlapping block vector by copying the
     *        front entries from their master process
     *
     * \todo use specialized send methods for improved
     *       performance. (i.e. only send the front entries to the
     *       peers.)
     */
    void syncFrontFromMaster()
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
            receiveFrontFromMaster_(peerRank);
        }
        
        // wait until we have send everything
        waitSendFinished_();
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

    /*!
     * \brief Set all remote entities to a given scalar value
     */
    void resetRemote(Scalar value = 0.0)
    {
        int numDomestic = overlap_->numDomestic();
        for (int i = overlap_->numLocal(); i < numDomestic; ++i) {
            (*this)[i] = value;
        };
    }

    void print() const 
    {
        for (int i = 0; i < this->size(); ++i) {
            std::cout << "row " << i << (overlap_->isLocal(i)?" ":"*") << ": " << (*this)[i] << "\n";
        };
    };

private:
    void createBuffers_()
    {
        typename PeerSet::const_iterator peerIt;
        typename PeerSet::const_iterator peerEndIt = overlap_->peerSet().end();
        
        // send all indices to the peers
        peerIt = overlap_->peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;

            const DomesticOverlapWithPeer &domesticOverlap = overlap_->domesticOverlapWithPeer(peerRank);
            int numEntries = domesticOverlap.size();
            indicesSendBuff_[peerRank] = std::shared_ptr<MpiBuffer<RowIndex> >(new MpiBuffer<RowIndex>(numEntries));
            valuesSendBuff_[peerRank] = std::shared_ptr<MpiBuffer<FieldVector> >(new MpiBuffer<FieldVector>(numEntries));
            
            // fill the indices buffer with global indices
            MpiBuffer<RowIndex> &indicesSendBuff = *indicesSendBuff_[peerRank];
            typename DomesticOverlapWithPeer::const_iterator domIt = domesticOverlap.begin();
            typename DomesticOverlapWithPeer::const_iterator domEndIt = domesticOverlap.end();
            for (int i = 0; domIt != domEndIt; ++domIt, ++i) {
                int rowIdx = *domIt;
                indicesSendBuff[i] = overlap_->domesticToGlobal(rowIdx);
            }

            // first, send the number of indices
            MPI_Bsend(&numEntries, // buff
                      1, // count
                      MPI_INT, // data type
                      peerRank, 
                      0, // tag
                      MPI_COMM_WORLD); // communicator
            // then, send the indices themselfs
            indicesSendBuff.send(peerRank);

            // change the indices of the send buffer to domestic ones
            domIt = domesticOverlap.begin();
            for (int i = 0; domIt != domEndIt; ++domIt, ++i) {
                int rowIdx = *domIt;
                indicesSendBuff[i] = rowIdx;
            }
        }

        // receive the indices from the peers
        peerIt = overlap_->peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;

            int numEntries;
            // first, receive the number of indices
            MPI_Recv(&numEntries, // buff
                     1, // count
                     MPI_INT, // data type
                     peerRank, 
                     0, // tag
                     MPI_COMM_WORLD, // communicator
                     MPI_STATUS_IGNORE);

            // then, create the MPI buffers
            indicesRecvBuff_[peerRank] = std::shared_ptr<MpiBuffer<RowIndex> >(new MpiBuffer<RowIndex>(numEntries));
            valuesRecvBuff_[peerRank] = std::shared_ptr<MpiBuffer<FieldVector> >(new MpiBuffer<FieldVector>(numEntries));
            MpiBuffer<RowIndex> &indicesRecvBuff = *indicesRecvBuff_[peerRank];

            // next, receive the actual indices
            indicesRecvBuff.receive(peerRank);

            // finally, translate the global indices to domestic ones
            for (int i = 0; i != numEntries; ++i) {
                int globalRowIdx = indicesRecvBuff[i];
                indicesRecvBuff[i] = overlap_->globalToDomestic(globalRowIdx);
            }
        }

        // wait for all send operations to complete
        peerIt = overlap_->peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            indicesSendBuff_[peerRank]->wait();
        }
    }

    void sendEntries_(int peerRank)
    {
        // copy the values into the send buffer
        const MpiBuffer<RowIndex> &indices = *indicesSendBuff_[peerRank];
        MpiBuffer<FieldVector> &values = *valuesSendBuff_[peerRank];
        for (int i = 0; i < indices.size(); ++ i)
            values[i] = (*this)[indices[i]];

        values.send(peerRank);
    }

    void waitSendFinished_()
    {
        typename PeerSet::const_iterator peerIt;
        typename PeerSet::const_iterator peerEndIt = overlap_->peerSet().end();

        // send all entries to all peers
        peerIt = overlap_->peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            valuesSendBuff_[peerRank]->wait();
        }
    }

    void receiveAddEntries_(int peerRank)
    {
        const MpiBuffer<RowIndex> &indices = *indicesRecvBuff_[peerRank];
        MpiBuffer<FieldVector> &values = *valuesRecvBuff_[peerRank];

        // receive the values from the peer
        values.receive(peerRank);

        // copy them into the block vector
        for (int i = 0; i < indices.size(); ++ i)
            (*this)[indices[i]] += values[i];
    }

    void receiveFrontFromMaster_(int peerRank)
    {
        const MpiBuffer<RowIndex> &indices = *indicesRecvBuff_[peerRank];
        MpiBuffer<FieldVector> &values = *valuesRecvBuff_[peerRank];

        // receive the values from the peer
        values.receive(peerRank);

        // copy them into the block vector
        for (int j = 0; j < indices.size(); ++j) {
            int domRowIdx = indices[j];
            if (overlap_->isFront(domRowIdx) && 
                overlap_->isMasterOf(peerRank, domRowIdx))
            {
                (*this)[domRowIdx] = values[j];
            }
        }
    }

    std::map<ProcessRank, std::shared_ptr<MpiBuffer<RowIndex> > > indicesSendBuff_;
    std::map<ProcessRank, std::shared_ptr<MpiBuffer<RowIndex> > > indicesRecvBuff_;

    std::map<ProcessRank, std::shared_ptr<MpiBuffer<FieldVector> > > valuesSendBuff_;
    std::map<ProcessRank, std::shared_ptr<MpiBuffer<FieldVector> > > valuesRecvBuff_;

    const Overlap *overlap_;
};

} // namespace Dumux

#endif
