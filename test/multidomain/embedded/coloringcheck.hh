// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief Check of the coloring of an embedded coupling manager for multithreaded assembly
 */
#ifndef DUMUX_TEST_MULTIDOMAIN_EMBEDDED_COLORINGCHECK_HH
#define DUMUX_TEST_MULTIDOMAIN_EMBEDDED_COLORINGCHECK_HH

#include <cstddef>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <utility>
#include <vector>

#include <dune/common/hybridutilities.hh>
#include <dune/common/indices.hh>

#include <dumux/discretization/method.hh>

namespace Dumux::Test {

/*!
 * \brief Count the pairs of elements with the same color whose assembly conflicts
 *
 * The assembly of an element deflects the solution at its dofs, at the dofs of its coupling stencils and,
 * in the bulk domain, at the dofs of its extended source stencil. Two elements conflict if one of them
 * deflects a dof the other one reads. An element reads its dofs, the dofs of its face neighbors if it is a
 * cell-centered element, and the dofs of its stencils. The check is independent of the coloring algorithm.
 */
template<class CouplingManager>
std::size_t countColoringConflicts(const CouplingManager& couplingManager)
{
    using MDTraits = typename CouplingManager::MultiDomainTraits;
    constexpr auto bulkIdx = Dune::index_constant<0>{};
    constexpr std::size_t numDomains = MDTraits::numSubDomains;

    std::size_t numConflicts = 0;
    using namespace Dune::Hybrid;
    forEach(integralRange(Dune::index_constant<numDomains>{}), [&](const auto i)
    {
        const auto& gg = couplingManager.problem(i).gridGeometry();
        const auto& grid = gg.gridView().grid();
        constexpr bool cellCentered = std::decay_t<decltype(gg)>::discMethod == DiscretizationMethods::cctpfa;

        // the dofs of all domains as one index
        const auto key = [&](std::size_t domain, std::size_t dof)
        { return dof*numDomains + domain; };

        const auto addDofs = [&](std::vector<std::size_t>& keys, const auto& element, bool deflected)
        {
            const auto eIdx = gg.elementMapper().index(element);
            if constexpr (cellCentered)
                keys.push_back(key(i, eIdx));
            else
                for (const auto dof : couplingManager.elementDofIndices(i, eIdx))
                    keys.push_back(key(i, dof));
            if constexpr (cellCentered)
                if (!deflected)
                    for (const auto& intersection : intersections(gg.gridView(), element))
                        if (intersection.neighbor())
                            keys.push_back(key(i, gg.elementMapper().index(intersection.outside())));
            forEach(integralRange(Dune::index_constant<numDomains>{}), [&](const auto j)
            {
                if constexpr (i != j)
                    for (const auto dof : couplingManager.couplingStencil(i, element, j))
                        keys.push_back(key(j, dof));
            });
            if constexpr (i == bulkIdx)
                for (const auto dof : couplingManager.extendedSourceStencil(eIdx))
                    keys.push_back(key(i, dof));
        };

        for (const auto& elementSet : couplingManager.elementSets(i))
        {
            std::unordered_map<std::size_t, std::pair<std::size_t, std::size_t>> deflectors; // dof -> (element, number)
            for (const auto& seed : elementSet)
            {
                const auto element = grid.entity(seed);
                const auto eIdx = gg.elementMapper().index(element);
                std::vector<std::size_t> keys;
                addDofs(keys, element, true);
                for (const auto k : keys)
                {
                    auto [it, inserted] = deflectors.try_emplace(k, std::make_pair(eIdx, std::size_t(1)));
                    if (!inserted && it->second.first != eIdx)
                    {
                        ++it->second.second;
                        ++numConflicts;
                    }
                }
            }

            for (const auto& seed : elementSet)
            {
                const auto element = grid.entity(seed);
                const auto eIdx = gg.elementMapper().index(element);
                std::vector<std::size_t> keys;
                addDofs(keys, element, false);
                for (const auto k : keys)
                {
                    const auto it = deflectors.find(k);
                    if (it != deflectors.end() && it->second.first != eIdx)
                        ++numConflicts;
                }
            }
        }
    });
    return numConflicts;
}

} // end namespace Dumux::Test

#endif
