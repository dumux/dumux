// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Parallel
 * \brief A dynamic sparse data exchange (DSDE) of byte messages between processes
 *
 * Implements the nonblocking-consensus (NBX) algorithm of Hoefler et al. (2010),
 * "Scalable Communication Protocols for Dynamic Sparse Data Exchange". Each process
 * knows only the messages it sends; the senders and message sizes are discovered
 * without any global (O(number of processes)) data structure or dense collective,
 * using nonblocking synchronous sends, message probing and a nonblocking barrier.
 */
#ifndef DUMUX_PARALLEL_NONBLOCKING_SPARSE_EXCHANGE_HH
#define DUMUX_PARALLEL_NONBLOCKING_SPARSE_EXCHANGE_HH

#include <vector>
#include <cstddef>
#include <cassert>

#if HAVE_MPI
#include <mpi.h>
#endif

namespace Dumux::Detail {

/*!
 * \ingroup Parallel
 * \brief Exchange variable-sized byte messages with a sparse set of peer processes
 * \param comm a Dune communicator (convertible to an MPI_Comm when built with MPI)
 * \param sendBufs one byte buffer per destination rank; an empty buffer means "no
 *                 message for that rank". A buffer addressed to the calling rank is ignored.
 * \param tag the MPI message tag to use; it must not clash with other communication
 *            in flight on the same communicator while this exchange runs.
 * \return the concatenation, in unspecified order, of all messages addressed to this rank
 *
 * \note Collective over \a comm. Each rank only needs to know its own destinations;
 *       who sends to it (and how much) is discovered on the fly, so there is no
 *       O(number of ranks) memory or all-to-all collective involved.
 */
template<class Communication>
std::vector<char> exchangeSparse(const Communication& comm,
                                 std::vector<std::vector<char>>& sendBufs,
                                 [[maybe_unused]] int tag = 7373)
{
    std::vector<char> recvBuf;

#if HAVE_MPI
    const int myRank = comm.rank();
    const int numProc = comm.size();
    if (numProc <= 1)
        return recvBuf;

    // Dune::Communication<MPI_Comm> exposes operator MPI_Comm() const
    const MPI_Comm mpiComm = comm;

    // 1. Post a nonblocking synchronous send for every non-empty message. An Issend
    //    completes locally only once the matching receive has been posted remotely,
    //    which is what lets the barrier below detect global completion.
    std::vector<MPI_Request> sendRequests;
    sendRequests.reserve(numProc);
    for (int to = 0; to < numProc; ++to)
    {
        if (to == myRank || sendBufs[to].empty())
            continue;
        sendRequests.emplace_back();
        MPI_Issend(sendBufs[to].data(), static_cast<int>(sendBufs[to].size()),
                   MPI_BYTE, to, tag, mpiComm, &sendRequests.back());
    }

    // Receive (via matched probe) every message that can currently be matched,
    // appending its bytes to the receive buffer.
    const auto receivePending = [&]
    {
        int flag = 0;
        MPI_Message message;
        MPI_Status status;
        MPI_Improbe(MPI_ANY_SOURCE, tag, mpiComm, &flag, &message, &status);
        while (flag)
        {
            int count = 0;
            MPI_Get_count(&status, MPI_BYTE, &count);
            const std::size_t offset = recvBuf.size();
            recvBuf.resize(offset + static_cast<std::size_t>(count));
            MPI_Mrecv(recvBuf.data() + offset, count, MPI_BYTE, &message, MPI_STATUS_IGNORE);
            MPI_Improbe(MPI_ANY_SOURCE, tag, mpiComm, &flag, &message, &status);
        }
    };

    // 2./3. Nonblocking consensus: keep draining incoming messages. Once all of our
    //       own sends have been matched we enter a nonblocking barrier; its completion
    //       means every rank's sends have been matched globally, so we are done.
    MPI_Request barrierRequest = MPI_REQUEST_NULL;
    bool barrierActive = false;
    bool done = false;
    while (!done)
    {
        receivePending();

        if (barrierActive)
        {
            int barrierDone = 0;
            MPI_Test(&barrierRequest, &barrierDone, MPI_STATUS_IGNORE);
            done = (barrierDone != 0);
        }
        else
        {
            int sendsDone = 0;
            MPI_Testall(static_cast<int>(sendRequests.size()),
                        sendRequests.data(), &sendsDone, MPI_STATUSES_IGNORE);
            if (sendsDone)
            {
                MPI_Ibarrier(mpiComm, &barrierRequest);
                barrierActive = true;
            }
        }
    }
#else
    assert(comm.size() == 1 && "Sparse data exchange between multiple ranks requires MPI");
    (void)sendBufs;
#endif

    return recvBuf;
}

} // end namespace Dumux::Detail

#endif
