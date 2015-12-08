// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
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

#warning This file is deprecated and will be removed after Dumux 2.9

#include <iostream>
#include <map>
#include <memory>
#include <vector>

#include <dune/istl/bvector.hh>

#include <dumux/parallel/mpibuffer.hh>
#include <dumux/common/valgrind.hh>

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
    typedef typename Overlap::DomesticOverlapWithPeer DomesticOverlapWithPeer;

    typedef typename ParentType::field_type Scalar;

public:
    OverlappingBlockVector(const Overlap &overlap)
        : ParentType(overlap.numDomestic())
        , overlap_(&overlap)
    {
        createBuffers_();
    }

    /*!
     * \brief Copy constructor.
     */
    OverlappingBlockVector(const OverlappingBlockVector &obv)
        : ParentType(obv)
        , frontMaster_(obv.frontMaster_)
        , numIndicesSendBuff_(obv.numIndicesSendBuff_)
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
        // assign the local rows
        int numLocal = overlap_->numLocal();
        for (int rowIdx = 0; rowIdx < numLocal; ++rowIdx) {
            if (overlap_->iAmMasterOf(rowIdx) ||
                overlap_->isShared(rowIdx))
            {
                (*this)[rowIdx] = nbv[rowIdx];
            }
            else
                (*this)[rowIdx] = 0;

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
    }

    /*!
     * \brief Syncronize the block vector by taking the arithmetic
     *        mean of all entries which are not on the front of some
     *        process
     *
     * \todo use specialized send methods for improved
     *       performance. (i.e. only send the front entries to the
     *       peers.)
     */
    void syncAverageFrontFromMaster()
    {
        // first, reset all of our front rows
        int numLocal = overlap_->numLocal();
        int numDomestic = overlap_->numDomestic();
        for (int i = numLocal; i < numDomestic; ++i) {
            if (overlap_->isFront(i))
                (*this)[i] = 0;
        }

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
            receiveAverageFrontFromMaster_(peerRank);
        }

        // divide each entry by the number of non-front processes
        for (int i = 0; i < numDomestic; ++i) {
            (*this)[i] /= overlap_->numNonFrontProcesses(i);
        }

        // wait until we have send everything
        waitSendFinished_();
    }

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
    }

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
    }

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
        }
    }

    /*!
     * \brief Set all remote entities to a given scalar value
     */
    void resetRemote(Scalar value = 0.0)
    {
        int numDomestic = overlap_->numDomestic();
        for (int i = overlap_->numLocal(); i < numDomestic; ++i) {
            (*this)[i] = value;
        }
    }

    void print() const
    {
        for (int i = 0; i < this->size(); ++i) {
            std::cout << "row " << i << (overlap_->isLocal(i)?" ":"*") << ": " << (*this)[i] << "\n";
        }
    }

private:
    void createBuffers_()
    {
#if HAVE_MPI
        // create array for the front indices
        int numDomestic = overlap_->numDomestic();
        frontMaster_ = std::make_shared<std::vector<ProcessRank> >(numDomestic, -1);

        typename PeerSet::const_iterator peerIt;
        typename PeerSet::const_iterator peerEndIt = overlap_->peerSet().end();

        // send all indices to the peers
        peerIt = overlap_->peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;

            const DomesticOverlapWithPeer &domesticOverlap = overlap_->domesticOverlapWithPeer(peerRank);
            int numEntries = domesticOverlap.size();
            numIndicesSendBuff_[peerRank] = std::make_shared<MpiBuffer<int> >(1);
            indicesSendBuff_[peerRank] = std::make_shared<MpiBuffer<RowIndex> >(numEntries);
            valuesSendBuff_[peerRank] = std::make_shared<MpiBuffer<FieldVector> >(numEntries);

            // fill the indices buffer with global indices
            MpiBuffer<RowIndex> &indicesSendBuff = *indicesSendBuff_[peerRank];
            typename DomesticOverlapWithPeer::const_iterator domIt = domesticOverlap.begin();
            typename DomesticOverlapWithPeer::const_iterator domEndIt = domesticOverlap.end();
            for (int i = 0; domIt != domEndIt; ++domIt, ++i) {
                int rowIdx = *domIt;
                indicesSendBuff[i] = overlap_->domesticToGlobal(rowIdx);
            }

            // first, send the number of indices
            (*numIndicesSendBuff_[peerRank])[0] = numEntries;
            numIndicesSendBuff_[peerRank]->send(peerRank);

            // then, send the indices themselfs
            indicesSendBuff.send(peerRank);
        }

        // receive the indices from the peers
        peerIt = overlap_->peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;

            // receive size of overlap to peer
            int numRows;
            MpiBuffer<int> numRowsRecvBuff(1);
            numRowsRecvBuff.receive(peerRank);
            numRows = numRowsRecvBuff[0];

            // then, create the MPI buffers
            indicesRecvBuff_[peerRank] = std::make_shared<MpiBuffer<RowIndex> >(numRows);
            valuesRecvBuff_[peerRank] = std::make_shared<MpiBuffer<FieldVector> >(numRows);
            MpiBuffer<RowIndex> &indicesRecvBuff = *indicesRecvBuff_[peerRank];

            // next, receive the actual indices
            indicesRecvBuff.receive(peerRank);

            // finally, translate the global indices to domestic ones
            for (int i = 0; i != numRows; ++i) {
                int globalRowIdx = indicesRecvBuff[i];
                int domRowIdx = overlap_->globalToDomestic(globalRowIdx);
                indicesRecvBuff[i] = domRowIdx;

                if (overlap_->isFront(domRowIdx) &&
                    overlap_->isMasterOf(peerRank, domRowIdx))
                {
                    (*frontMaster_)[domRowIdx] = peerRank;
                }
            }
        }

        // wait for all send operations to complete
        peerIt = overlap_->peerSet().begin();
        for (; peerIt != peerEndIt; ++peerIt) {
            int peerRank = *peerIt;
            numIndicesSendBuff_[peerRank]->wait();
            indicesSendBuff_[peerRank]->wait();

            // convert the global indices of the send buffer to
            // domestic ones
            MpiBuffer<RowIndex> &indicesSendBuff = *indicesSendBuff_[peerRank];
            for (int i = 0; i < indicesSendBuff.size(); ++i) {
                indicesSendBuff[i] = overlap_->globalToDomestic(indicesSendBuff[i]);
            }
        }
#endif // HAVE_MPI
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
        for (int i = 0; i < indices.size(); ++ i) {
            (*this)[indices[i]] += values[i];
        }
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
            if ((*frontMaster_)[domRowIdx] != peerRank)
                continue;

            (*this)[domRowIdx] = values[j];
        }
    }

    void receiveAverageFrontFromMaster_(int peerRank)
    {
        const MpiBuffer<RowIndex> &indices = *indicesRecvBuff_[peerRank];
        MpiBuffer<FieldVector> &values = *valuesRecvBuff_[peerRank];

        // receive the values from the peer
        values.receive(peerRank);

        // copy them into the block vector
        for (int j = 0; j < indices.size(); ++j) {
            int domRowIdx = indices[j];
            (*this)[domRowIdx] += values[j];
        }
    }

    std::shared_ptr<std::vector<ProcessRank> > frontMaster_;

    std::map<ProcessRank, std::shared_ptr<MpiBuffer<RowIndex> > > numIndicesSendBuff_;
    std::map<ProcessRank, std::shared_ptr<MpiBuffer<RowIndex> > > indicesSendBuff_;
    std::map<ProcessRank, std::shared_ptr<MpiBuffer<RowIndex> > > indicesRecvBuff_;

    std::map<ProcessRank, std::shared_ptr<MpiBuffer<FieldVector> > > valuesSendBuff_;
    std::map<ProcessRank, std::shared_ptr<MpiBuffer<FieldVector> > > valuesRecvBuff_;

    const Overlap *overlap_;
};

} // namespace Dumux

#endif
