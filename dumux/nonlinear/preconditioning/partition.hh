// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Nonlinear
 * \brief Overlapping decomposition of a grid into subdomains for nonlinear domain decomposition
 */
#ifndef DUMUX_NONLINEAR_PRECONDITIONING_PARTITION_HH
#define DUMUX_NONLINEAR_PRECONDITIONING_PARTITION_HH

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>

namespace Dumux::NonlinearPreconditioning {

using SubdomainIndex = std::size_t;
inline constexpr SubdomainIndex noSubdomain = std::numeric_limits<SubdomainIndex>::max();

/*!
 * \ingroup Nonlinear
 * \brief An overlapping decomposition of a set of degrees of freedom into subdomains
 *
 * The decomposition is built from an influence map: `influence[i]` lists the degrees of freedom
 * whose residual depends on the degree of freedom `i`. This is exactly the structure
 * `CCSimpleConnectivityMap` provides, and exactly the structure that determines which matrix rows
 * an element contributes to when it is visited by an assembly loop.
 *
 * Three index sets are exported per subdomain:
 *  - the core, a disjoint cover of all degrees of freedom, defining the partition of unity;
 *  - the subdomain \f$ N_i \f$, the core grown by `overlapWidth` rings, whose rows and columns form
 *    the square local block used by the local nonlinear solve;
 *  - the fringe \f$ \Gamma_i \f$, those degrees of freedom outside \f$ N_i \f$ that still contribute
 *    to a row in \f$ N_i \f$. Visiting \f$ N_i \cup \Gamma_i \f$ yields the rectangular block with
 *    global columns that the exact RASPEN Jacobian requires.
 *
 * \note For cell-centred schemes an element index and its degree of freedom index coincide, so the
 *       same sets serve as assembly loop ranges and as matrix index sets. Vertex-centred schemes
 *       need a separate element-to-dof mapping and are not supported here.
 */
class Partition
{
public:
    using Index = std::size_t;
    static constexpr Index noIndex = std::numeric_limits<Index>::max();

    /*!
     * \brief Build a decomposition from an influence map and a disjoint core assignment
     * \param influence for each dof, the dofs whose residual depends on it
     * \param coreAssignment for each dof, the subdomain owning it in its core
     * \param overlapWidth number of rings the core is grown by
     */
    Partition(const std::vector<std::vector<Index>>& influence,
              std::vector<SubdomainIndex> coreAssignment,
              std::size_t overlapWidth)
    : coreAssignment_(std::move(coreAssignment))
    , overlapWidth_(overlapWidth)
    {
        if (influence.size() != coreAssignment_.size())
            DUNE_THROW(Dune::InvalidStateException,
                       "Influence map has " << influence.size() << " entries but the core assignment has "
                       << coreAssignment_.size());

        if (std::any_of(coreAssignment_.begin(), coreAssignment_.end(),
                        [] (auto s) { return s == noSubdomain; }))
            DUNE_THROW(Dune::InvalidStateException, "Core assignment does not cover all degrees of freedom");

        numSubdomains_ = coreAssignment_.empty()
            ? 0 : 1 + *std::max_element(coreAssignment_.begin(), coreAssignment_.end());

        buildCores_();
        buildSubdomains_(symmetrized_(influence));
        buildFringes_(influence);
        buildLocalIndices_();
        buildWeights_();
    }

    std::size_t numSubdomains() const { return numSubdomains_; }
    std::size_t overlapWidth() const { return overlapWidth_; }
    std::size_t numDofs() const { return coreAssignment_.size(); }

    SubdomainIndex coreSubdomain(Index dof) const { return coreAssignment_[dof]; }

    //! the disjoint cover; defines the partition of unity
    const std::vector<Index>& core(SubdomainIndex i) const { return core_[i]; }

    //! \f$ N_i \f$, the core grown by overlapWidth() rings
    const std::vector<Index>& dofs(SubdomainIndex i) const { return dofs_[i]; }

    //! \f$ \Gamma_i \f$, dofs outside \f$ N_i \f$ that contribute to a row in \f$ N_i \f$
    const std::vector<Index>& fringe(SubdomainIndex i) const { return fringe_[i]; }

    //! the assembly loop range yielding the rectangular block with global columns
    std::vector<Index> dofsWithFringe(SubdomainIndex i) const
    {
        std::vector<Index> result;
        result.reserve(dofs_[i].size() + fringe_[i].size());
        std::set_union(dofs_[i].begin(), dofs_[i].end(),
                       fringe_[i].begin(), fringe_[i].end(),
                       std::back_inserter(result));
        return result;
    }

    //! position of a global dof within dofs(i), or noIndex if it is not in \f$ N_i \f$
    Index localIndex(SubdomainIndex i, Index globalDof) const { return localIndex_[i][globalDof]; }

    //! \f$ D_i \f$, aligned with dofs(i); the indicator of the core
    const std::vector<double>& weights(SubdomainIndex i) const { return weights_[i]; }

private:
    //! overlap grows along a symmetric neighbourhood, whereas the fringe follows the write direction
    std::vector<std::vector<Index>> symmetrized_(const std::vector<std::vector<Index>>& influence) const
    {
        std::vector<std::vector<Index>> neighbours(influence.size());
        for (Index i = 0; i < influence.size(); ++i)
            for (const auto j : influence[i])
            {
                neighbours[i].push_back(j);
                neighbours[j].push_back(i);
            }

        for (auto& n : neighbours)
        {
            std::sort(n.begin(), n.end());
            n.erase(std::unique(n.begin(), n.end()), n.end());
        }

        return neighbours;
    }

    void buildCores_()
    {
        core_.assign(numSubdomains_, {});
        for (Index dof = 0; dof < coreAssignment_.size(); ++dof)
            core_[coreAssignment_[dof]].push_back(dof);
    }

    void buildSubdomains_(const std::vector<std::vector<Index>>& neighbours)
    {
        dofs_.assign(numSubdomains_, {});
        std::vector<bool> inSubdomain(coreAssignment_.size(), false);

        for (SubdomainIndex i = 0; i < numSubdomains_; ++i)
        {
            auto& set = dofs_[i];
            set = core_[i];
            for (const auto dof : set)
                inSubdomain[dof] = true;

            std::size_t ringBegin = 0;
            for (std::size_t ring = 0; ring < overlapWidth_; ++ring)
            {
                const auto ringEnd = set.size();
                for (std::size_t k = ringBegin; k < ringEnd; ++k)
                    for (const auto j : neighbours[set[k]])
                        if (!inSubdomain[j])
                        {
                            inSubdomain[j] = true;
                            set.push_back(j);
                        }
                ringBegin = ringEnd;
            }

            std::sort(set.begin(), set.end());
            for (const auto dof : set)
                inSubdomain[dof] = false;
        }
    }

    void buildFringes_(const std::vector<std::vector<Index>>& influence)
    {
        fringe_.assign(numSubdomains_, {});
        std::vector<bool> inSubdomain(coreAssignment_.size(), false);

        for (SubdomainIndex i = 0; i < numSubdomains_; ++i)
        {
            for (const auto dof : dofs_[i])
                inSubdomain[dof] = true;

            for (Index dof = 0; dof < influence.size(); ++dof)
                if (!inSubdomain[dof]
                    && std::any_of(influence[dof].begin(), influence[dof].end(),
                                   [&] (auto j) { return inSubdomain[j]; }))
                    fringe_[i].push_back(dof);

            for (const auto dof : dofs_[i])
                inSubdomain[dof] = false;
        }
    }

    void buildLocalIndices_()
    {
        localIndex_.assign(numSubdomains_, std::vector<Index>(coreAssignment_.size(), noIndex));
        for (SubdomainIndex i = 0; i < numSubdomains_; ++i)
            for (Index k = 0; k < dofs_[i].size(); ++k)
                localIndex_[i][dofs_[i][k]] = k;
    }

    void buildWeights_()
    {
        weights_.assign(numSubdomains_, {});
        for (SubdomainIndex i = 0; i < numSubdomains_; ++i)
        {
            weights_[i].reserve(dofs_[i].size());
            for (const auto dof : dofs_[i])
                weights_[i].push_back(coreAssignment_[dof] == i ? 1.0 : 0.0);
        }
    }

    std::vector<SubdomainIndex> coreAssignment_;
    std::size_t overlapWidth_;
    std::size_t numSubdomains_;

    std::vector<std::vector<Index>> core_;
    std::vector<std::vector<Index>> dofs_;
    std::vector<std::vector<Index>> fringe_;
    std::vector<std::vector<Index>> localIndex_;
    std::vector<std::vector<double>> weights_;
};

/*!
 * \ingroup Nonlinear
 * \brief Structural checks on a Partition
 */
struct PartitionValidator
{
    struct Report
    {
        bool covering = false;          //!< every dof lies in exactly one core
        bool partitionOfUnity = false;  //!< \f$ \sum_i R_i^T D_i R_i = I \f$
        bool coresConnected = false;    //!< each core is connected within its subdomain
        bool overlapWithinWidth = false;//!< every overlap dof is within overlapWidth() hops of the core

        bool valid() const { return covering && partitionOfUnity && coresConnected && overlapWithinWidth; }
    };

    static Report validate(const Partition& p, const std::vector<std::vector<Partition::Index>>& influence)
    {
        using Index = Partition::Index;
        Report report;

        std::size_t coreSize = 0;
        for (SubdomainIndex i = 0; i < p.numSubdomains(); ++i)
            coreSize += p.core(i).size();
        report.covering = (coreSize == p.numDofs());

        std::vector<double> weightSum(p.numDofs(), 0.0);
        for (SubdomainIndex i = 0; i < p.numSubdomains(); ++i)
            for (Index k = 0; k < p.dofs(i).size(); ++k)
                weightSum[p.dofs(i)[k]] += p.weights(i)[k];
        report.partitionOfUnity = std::all_of(weightSum.begin(), weightSum.end(),
                                              [] (auto w) { return std::abs(w - 1.0) < 1e-14; });

        std::vector<std::vector<Index>> neighbours(influence.size());
        for (Index i = 0; i < influence.size(); ++i)
            for (const auto j : influence[i])
            {
                neighbours[i].push_back(j);
                neighbours[j].push_back(i);
            }

        report.coresConnected = true;
        report.overlapWithinWidth = true;
        for (SubdomainIndex i = 0; i < p.numSubdomains(); ++i)
        {
            if (p.core(i).empty())
                continue;

            std::vector<bool> inCore(p.numDofs(), false);
            for (const auto dof : p.core(i))
                inCore[dof] = true;

            std::vector<bool> seen(p.numDofs(), false);
            std::vector<Index> queue{p.core(i).front()};
            seen[queue.front()] = true;
            for (std::size_t head = 0; head < queue.size(); ++head)
                for (const auto j : neighbours[queue[head]])
                    if (inCore[j] && !seen[j])
                    {
                        seen[j] = true;
                        queue.push_back(j);
                    }
            if (queue.size() != p.core(i).size())
                report.coresConnected = false;

            std::vector<std::size_t> distance(p.numDofs(), std::numeric_limits<std::size_t>::max());
            std::vector<Index> ring = p.core(i);
            for (const auto dof : ring)
                distance[dof] = 0;
            for (std::size_t d = 1; d <= p.overlapWidth(); ++d)
            {
                std::vector<Index> next;
                for (const auto dof : ring)
                    for (const auto j : neighbours[dof])
                        if (distance[j] == std::numeric_limits<std::size_t>::max())
                        {
                            distance[j] = d;
                            next.push_back(j);
                        }
                ring = std::move(next);
            }
            for (const auto dof : p.dofs(i))
                if (distance[dof] == std::numeric_limits<std::size_t>::max())
                    report.overlapWithinWidth = false;
        }

        return report;
    }
};

} // end namespace Dumux::NonlinearPreconditioning

#endif
