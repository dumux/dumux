// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Parallel
 * \brief Test the byte-buffer (de)serialization helpers and the nonblocking-consensus
 *        sparse data exchange.
 */
#include <config.h>

#include <set>
#include <vector>
#include <cstddef>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/initialize.hh>
#include <dumux/parallel/packing.hh>
#include <dumux/parallel/nonblockingsparseexchange.hh>

namespace {

//! Pack values of different trivially-copyable types and read them back
void testPacking()
{
    using namespace Dumux::Detail;

    const int i = -42;
    const std::size_t s = 1234567u;
    const double d = 3.14159265358979;
    const char ch = 'x';

    std::vector<char> buf;
    packValue(buf, i);
    packValue(buf, s);
    packValue(buf, d);
    packValue(buf, ch);

    if (buf.size() != sizeof(i) + sizeof(s) + sizeof(d) + sizeof(ch))
        DUNE_THROW(Dune::Exception, "packValue produced an unexpected buffer size");

    const char* cursor = buf.data();
    if (unpackValue<int>(cursor) != i)
        DUNE_THROW(Dune::Exception, "unpackValue<int> round trip failed");
    if (unpackValue<std::size_t>(cursor) != s)
        DUNE_THROW(Dune::Exception, "unpackValue<std::size_t> round trip failed");
    if (Dune::FloatCmp::ne(unpackValue<double>(cursor), d))
        DUNE_THROW(Dune::Exception, "unpackValue<double> round trip failed");
    if (unpackValue<char>(cursor) != ch)
        DUNE_THROW(Dune::Exception, "unpackValue<char> round trip failed");

    // the cursor must have advanced to exactly the end of the buffer
    if (cursor != buf.data() + buf.size())
        DUNE_THROW(Dune::Exception, "unpacking did not advance the cursor to the buffer end");
}

// deterministic, sparse send pattern: rank r sends to r+1 and r+2 (modulo size,
// skipping self), so for size <= 1 nobody sends and for size == 2 each rank sends once.
std::set<int> destinationsOf(int sender, int size)
{
    std::set<int> dests;
    for (const int offset : {1, 2})
    {
        const int to = (sender + offset) % size;
        if (to != sender)
            dests.insert(to);
    }
    return dests;
}

//! a sender-dependent, variable payload length to also exercise variable-size messages
int payloadLength(int sender)
{ return (sender % 3) + 1; }

//! Exchange a sparse set of self-describing messages and verify each rank gets exactly its own
template<class Communication>
void testExchange(const Communication& comm)
{
    using namespace Dumux::Detail;

    const int myRank = comm.rank();
    const int numProc = comm.size();

    // pack one message per destination: [sender][length][length x sender]
    std::vector<std::vector<char>> sendBufs(numProc);
    for (const int to : destinationsOf(myRank, numProc))
    {
        auto& buf = sendBufs[to];
        const int length = payloadLength(myRank);
        packValue(buf, myRank);
        packValue(buf, length);
        for (int k = 0; k < length; ++k)
            packValue(buf, myRank);
    }

    const auto recvBuf = exchangeSparse(comm, sendBufs);

    // the senders we expect to hear from
    std::set<int> expectedSenders;
    for (int s = 0; s < numProc; ++s)
        if (destinationsOf(s, numProc).count(myRank))
            expectedSenders.insert(s);

    // parse the concatenated messages and check their integrity
    std::set<int> gotSenders;
    const char* cursor = recvBuf.data();
    const char* const end = recvBuf.data() + recvBuf.size();
    while (cursor < end)
    {
        const int sender = unpackValue<int>(cursor);
        const int length = unpackValue<int>(cursor);
        if (length != payloadLength(sender))
            DUNE_THROW(Dune::Exception, "received message has unexpected payload length");
        for (int k = 0; k < length; ++k)
            if (unpackValue<int>(cursor) != sender)
                DUNE_THROW(Dune::Exception, "received message payload is corrupted");
        if (!gotSenders.insert(sender).second)
            DUNE_THROW(Dune::Exception, "received more than one message from the same sender");
    }

    if (cursor != end)
        DUNE_THROW(Dune::Exception, "received buffer has trailing/truncated bytes");
    if (gotSenders != expectedSenders)
        DUNE_THROW(Dune::Exception, "the set of senders does not match the expected sparse pattern");
}

} // end anonymous namespace

int main(int argc, char** argv)
{
    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    testPacking();
    testExchange(mpiHelper.getCommunication());

    if (mpiHelper.rank() == 0)
        std::cout << "All sparse-exchange tests passed on " << mpiHelper.size() << " rank(s)" << std::endl;

    return 0;
}
