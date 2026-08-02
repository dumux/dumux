// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>

#include <algorithm>
#include <iostream>
#include <vector>

#include <dumux/nonlinear/preconditioning/partition.hh>

namespace {

using Dumux::NonlinearPreconditioning::Partition;
using Dumux::NonlinearPreconditioning::SubdomainIndex;
using Index = Partition::Index;

int failures = 0;

void check(bool condition, const std::string& what)
{
    if (!condition)
    {
        std::cerr << "FAIL: " << what << std::endl;
        ++failures;
    }
}

void checkEqual(const std::vector<Index>& got, const std::vector<Index>& expected, const std::string& what)
{
    if (got != expected)
    {
        std::cerr << "FAIL: " << what << "\n  got     :";
        for (const auto i : got) std::cerr << " " << i;
        std::cerr << "\n  expected:";
        for (const auto i : expected) std::cerr << " " << i;
        std::cerr << std::endl;
        ++failures;
    }
}

//! influence map of a 1d chain of n cells with a two-point stencil
std::vector<std::vector<Index>> chainInfluence(std::size_t n)
{
    std::vector<std::vector<Index>> influence(n);
    for (Index i = 0; i + 1 < n; ++i)
    {
        influence[i].push_back(i + 1);
        influence[i + 1].push_back(i);
    }
    return influence;
}

//! influence map of an nx times ny structured grid with a two-point stencil
std::vector<std::vector<Index>> structuredInfluence(std::size_t nx, std::size_t ny)
{
    std::vector<std::vector<Index>> influence(nx*ny);
    const auto idx = [nx] (std::size_t i, std::size_t j) { return j*nx + i; };
    for (std::size_t j = 0; j < ny; ++j)
        for (std::size_t i = 0; i < nx; ++i)
        {
            if (i + 1 < nx) { influence[idx(i,j)].push_back(idx(i+1,j)); influence[idx(i+1,j)].push_back(idx(i,j)); }
            if (j + 1 < ny) { influence[idx(i,j)].push_back(idx(i,j+1)); influence[idx(i,j+1)].push_back(idx(i,j)); }
        }
    return influence;
}

void testChain()
{
    const auto influence = chainInfluence(10);
    const std::vector<SubdomainIndex> cores{0,0,0,0,0,1,1,1,1,1};

    {
        const Partition p(influence, cores, 0);
        checkEqual(p.dofs(0), {0,1,2,3,4}, "chain, overlap 0: N_0");
        checkEqual(p.fringe(0), {5}, "chain, overlap 0: fringe_0");
        checkEqual(p.fringe(1), {4}, "chain, overlap 0: fringe_1");
    }

    {
        const Partition p(influence, cores, 1);
        check(p.numSubdomains() == 2, "chain: two subdomains");
        checkEqual(p.dofs(0), {0,1,2,3,4,5}, "chain, overlap 1: N_0");
        checkEqual(p.dofs(1), {4,5,6,7,8,9}, "chain, overlap 1: N_1");
        checkEqual(p.fringe(0), {6}, "chain, overlap 1: fringe_0");
        checkEqual(p.fringe(1), {3}, "chain, overlap 1: fringe_1");
        checkEqual(p.dofsWithFringe(0), {0,1,2,3,4,5,6}, "chain, overlap 1: N_0 union fringe_0");

        // the core indicator is the partition of unity: cell 5 is overlap for subdomain 0, core for 1
        check(p.weights(0)[p.localIndex(0,5)] == 0.0, "chain: overlap dof has zero weight in the non-owning subdomain");
        check(p.weights(1)[p.localIndex(1,5)] == 1.0, "chain: core dof has unit weight in the owning subdomain");
        check(p.localIndex(0, 9) == Partition::noIndex, "chain: dof outside N_0 has no local index");

        const auto report = Dumux::NonlinearPreconditioning::PartitionValidator::validate(p, influence);
        check(report.valid(), "chain, overlap 1: partition is structurally valid");
    }

    {
        const Partition p(influence, cores, 2);
        checkEqual(p.dofs(0), {0,1,2,3,4,5,6}, "chain, overlap 2: N_0");
        checkEqual(p.fringe(0), {7}, "chain, overlap 2: fringe_0");
    }
}

void testStructured()
{
    // a 5x5 grid split into a left (i<3) and a right (i>=3) subdomain
    const std::size_t nx = 5, ny = 5;
    const auto influence = structuredInfluence(nx, ny);
    std::vector<SubdomainIndex> cores(nx*ny);
    for (std::size_t j = 0; j < ny; ++j)
        for (std::size_t i = 0; i < nx; ++i)
            cores[j*nx + i] = (i < 3) ? 0 : 1;

    const Partition p(influence, cores, 1);

    std::vector<Index> expected;
    for (std::size_t j = 0; j < ny; ++j)
        for (std::size_t i = 0; i < 4; ++i)   // core columns 0-2 grown by one column
            expected.push_back(j*nx + i);
    checkEqual(p.dofs(0), expected, "5x5, overlap 1: N_0 is the first four columns");

    expected.clear();
    for (std::size_t j = 0; j < ny; ++j)
        expected.push_back(j*nx + 4);
    checkEqual(p.fringe(0), expected, "5x5, overlap 1: fringe_0 is the last column");

    const auto report = Dumux::NonlinearPreconditioning::PartitionValidator::validate(p, influence);
    check(report.valid(), "5x5: partition is structurally valid");
}

//! the identity every RASPEN weight must satisfy, on an irregular partition
void testPartitionOfUnity()
{
    const auto influence = structuredInfluence(7, 6);
    std::vector<SubdomainIndex> cores(7*6);
    for (Index k = 0; k < cores.size(); ++k)
        cores[k] = (k % 5 == 0) ? 2 : (k < 20 ? 0 : 1);

    for (std::size_t overlap : {0ul, 1ul, 2ul, 3ul})
    {
        const Partition p(influence, cores, overlap);
        std::vector<double> sum(p.numDofs(), 0.0);
        for (SubdomainIndex i = 0; i < p.numSubdomains(); ++i)
            for (Index k = 0; k < p.dofs(i).size(); ++k)
                sum[p.dofs(i)[k]] += p.weights(i)[k];

        check(std::all_of(sum.begin(), sum.end(), [] (auto w) { return std::abs(w - 1.0) < 1e-14; }),
              "partition of unity holds at overlap " + std::to_string(overlap));
    }
}

//! every dof writing into a row of N_i must be in N_i or in the fringe, or the local block is incomplete
void testFringeCompleteness()
{
    const auto influence = structuredInfluence(6, 6);
    std::vector<SubdomainIndex> cores(36);
    for (std::size_t j = 0; j < 6; ++j)
        for (std::size_t i = 0; i < 6; ++i)
            cores[j*6 + i] = (i < 3 ? 0u : 1u) + (j < 3 ? 0u : 2u);

    for (std::size_t overlap : {0ul, 1ul, 2ul})
    {
        const Partition p(influence, cores, overlap);
        for (SubdomainIndex i = 0; i < p.numSubdomains(); ++i)
        {
            const auto range = p.dofsWithFringe(i);
            std::vector<bool> inRange(p.numDofs(), false);
            for (const auto dof : range)
                inRange[dof] = true;

            std::vector<bool> inSubdomain(p.numDofs(), false);
            for (const auto dof : p.dofs(i))
                inSubdomain[dof] = true;

            for (Index dof = 0; dof < influence.size(); ++dof)
            {
                const bool writesIntoSubdomain
                    = inSubdomain[dof]
                      || std::any_of(influence[dof].begin(), influence[dof].end(),
                                     [&] (auto j) { return inSubdomain[j]; });
                if (writesIntoSubdomain)
                    check(inRange[dof], "fringe completeness at overlap " + std::to_string(overlap)
                                        + ", subdomain " + std::to_string(i)
                                        + ": dof " + std::to_string(dof) + " must be visited");
            }
        }
    }
}

} // end anonymous namespace

int main()
{
    testChain();
    testStructured();
    testPartitionOfUnity();
    testFringeCompleteness();

    if (failures > 0)
    {
        std::cerr << failures << " check(s) failed" << std::endl;
        return 1;
    }

    std::cout << "All partition checks passed" << std::endl;
    return 0;
}
