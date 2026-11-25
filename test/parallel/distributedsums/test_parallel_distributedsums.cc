// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test the sums whose terms are computed on several processes
 *
 * Every rank p creates sums (p, k). Rank r holds a term of sum (p, k) if r == p or (r + p + k) is even,
 * and needs the sum if r == p or r == p + 1 (modulo the number of ranks). Ranks insert the sums in
 * different orders. Every rank that needs a sum compares its own term plus the received terms with the
 * sum computed from the formula for the terms, for two different sets of terms.
 */
#include <config.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <span>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/initialize.hh>
#include <dumux/parallel/distributedsums.hh>

namespace {

constexpr int sumsPerRank = 5;

std::size_t numComponents(int k)
{ return 1 + k % 3; }

bool holdsTerm(int r, int p, int k)
{ return r == p || (r + p + k) % 2 == 0; }

bool needsSum(int r, int p, int numRanks)
{ return r == p || r == (p + 1) % numRanks; }

double term(int r, int p, int k, std::size_t c, int round)
{ return (r + 1)*(p + 2)*(k + 3)*(c + 1) + 0.25*round*r; }

} // end anonymous namespace

int main(int argc, char** argv)
{
    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();
    const auto comm = mpiHelper.getCommunication();
    const int rank = comm.rank();
    const int numRanks = comm.size();

    Dumux::DistributedSums<double> sums(comm);

    std::vector<std::pair<int, int>> order;
    for (int p = 0; p < numRanks; ++p)
        for (int k = 0; k < sumsPerRank; ++k)
            order.emplace_back(p, k);
    if (rank % 2 == 1)
        std::reverse(order.begin(), order.end());

    for (const auto& [p, k] : order)
    {
        const bool holds = holdsTerm(rank, p, k);
        const bool needs = needsSum(rank, p, numRanks);
        if (!holds && !needs)
            continue;

        const auto i = sums.insert({p, static_cast<std::size_t>(k)}, numComponents(k));
        for (int r = 0; r < numRanks; ++r)
        {
            if (holds && needsSum(r, p, numRanks))
                sums.addDestination(i, r);
            if (needs && holdsTerm(r, p, k))
                sums.addSource(i, r);
        }
    }

    if (sums.index({rank, sumsPerRank}) != Dumux::DistributedSums<double>::invalidIndex)
        DUNE_THROW(Dune::Exception, "A sum that was never added has an index");

    sums.finalize();

    int failures = 0;
    for (int round = 0; round < 2; ++round)
    {
        sums.exchange([&](std::size_t i, std::span<double> values)
        {
            const auto [p, k] = sums.key(i);
            for (std::size_t c = 0; c < values.size(); ++c)
                values[c] = term(rank, p, k, c, round);
        });

        for (std::size_t i = 0; i < sums.size(); ++i)
        {
            const auto [p, k] = sums.key(i);
            if (!needsSum(rank, p, numRanks))
                continue;

            const auto received = sums.received(i);
            for (std::size_t c = 0; c < numComponents(k); ++c)
            {
                double expected = 0.0;
                for (int r = 0; r < numRanks; ++r)
                    if (holdsTerm(r, p, k))
                        expected += term(r, p, k, c, round);

                const double own = holdsTerm(rank, p, k) ? term(rank, p, k, c, round) : 0.0;
                if (std::abs(own + received[c] - expected) > 1e-12*std::abs(expected))
                {
                    std::cerr << "rank " << rank << " round " << round << ": sum (" << p << ", " << k
                              << ") component " << c << " is " << own + received[c]
                              << " instead of " << expected << std::endl;
                    ++failures;
                }
            }
        }
    }

    failures = comm.sum(failures);
    if (failures > 0)
        DUNE_THROW(Dune::Exception, failures << " wrong sums");

    if (rank == 0)
        std::cout << "distributed sums: all sums correct on " << numRanks << " rank(s)" << std::endl;

    return 0;
}
