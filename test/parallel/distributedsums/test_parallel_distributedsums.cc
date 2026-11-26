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
 * sum computed from the formula for the terms, for two different sets of terms, and for a third set
 * received into separate values, which must not change the received sums. Ranks also send lists of
 * different lengths, some empty, along the sums; a rank that needs a sum receives the lists of the other
 * ranks holding a term, in the order of their ranks.
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

struct ListItem
{
    int rank;
    int p;
    int k;
    int j;

    bool operator==(const ListItem&) const = default;
};

//! the list of rank r for sum (p, k)
std::vector<ListItem> list(int r, int p, int k)
{
    std::vector<ListItem> items;
    for (int j = 0; j < (r + k) % 3; ++j)
        items.push_back({r, p, k, j});
    return items;
}

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

    std::vector<double> values;
    sums.exchange([&](std::size_t i, std::span<double> v)
    {
        const auto [p, k] = sums.key(i);
        for (std::size_t c = 0; c < v.size(); ++c)
            v[c] = term(rank, p, k, c, 2);
    }, values);

    for (std::size_t i = 0; i < sums.size(); ++i)
    {
        const auto [p, k] = sums.key(i);
        bool othersHoldTerms = false;
        for (int r = 0; r < numRanks; ++r)
            othersHoldTerms = othersHoldTerms || (r != rank && holdsTerm(r, p, k));
        if (sums.receives(i) != (needsSum(rank, p, numRanks) && othersHoldTerms))
        {
            std::cerr << "rank " << rank << ": receives() of sum (" << p << ", " << k << ") is wrong" << std::endl;
            ++failures;
        }

        if (!needsSum(rank, p, numRanks))
            continue;

        for (std::size_t c = 0; c < numComponents(k); ++c)
        {
            double expectedValues = 0.0, expectedReceived = 0.0;
            for (int r = 0; r < numRanks; ++r)
            {
                if (r != rank && holdsTerm(r, p, k))
                {
                    expectedValues += term(r, p, k, c, 2);
                    expectedReceived += term(r, p, k, c, 1);
                }
            }
            if (std::abs(sums.received(i, values)[c] - expectedValues) > 1e-12*std::abs(expectedValues)
                || std::abs(sums.received(i)[c] - expectedReceived) > 1e-12*std::abs(expectedReceived))
            {
                std::cerr << "rank " << rank << ": sum (" << p << ", " << k << ") component " << c
                          << " received into separate values is wrong or changed the received sums" << std::endl;
                ++failures;
            }
        }
    }

    const auto lists = sums.exchangeLists<ListItem>([&](std::size_t i)
    {
        const auto [p, k] = sums.key(i);
        return list(rank, p, k);
    });
    for (std::size_t i = 0; i < sums.size(); ++i)
    {
        const auto [p, k] = sums.key(i);
        std::vector<ListItem> expected;
        if (needsSum(rank, p, numRanks))
            for (int r = 0; r < numRanks; ++r)
                if (r != rank && holdsTerm(r, p, k))
                    for (const auto& item : list(r, p, k))
                        expected.push_back(item);
        if (lists[i] != expected)
        {
            std::cerr << "rank " << rank << ": the lists of sum (" << p << ", " << k << ") are wrong" << std::endl;
            ++failures;
        }
    }

    failures = comm.sum(failures);
    if (failures > 0)
        DUNE_THROW(Dune::Exception, failures << " wrong sums");

    if (rank == 0)
        std::cout << "distributed sums: all sums correct on " << numRanks << " rank(s)" << std::endl;

    return 0;
}
