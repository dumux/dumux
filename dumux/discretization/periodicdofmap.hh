// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Helpers for maps between degrees of freedom identified across periodic boundaries
 */
#ifndef DUMUX_DISCRETIZATION_PERIODIC_DOF_MAP_HH
#define DUMUX_DISCRETIZATION_PERIODIC_DOF_MAP_HH

#include <algorithm>
#include <cstddef>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Dumux::Detail {

/*!
 * \ingroup Discretization
 * \brief Register a degree of freedom as periodically mapped to another one
 */
template<class Index>
void addPeriodicallyMappedDof(std::unordered_map<Index, std::vector<Index>>& periodicDofMap,
                              std::type_identity_t<Index> dofIdx, std::type_identity_t<Index> periodicDofIdx)
{
    auto& periodicDofs = periodicDofMap[dofIdx];
    if (std::ranges::find(periodicDofs, periodicDofIdx) == periodicDofs.end())
        periodicDofs.push_back(periodicDofIdx);
}

/*!
 * \ingroup Discretization
 * \brief Map every degree of freedom to all others it is identified with
 *
 * The map of degrees of freedom to their partners across a single periodic boundary
 * is extended to the transitive closure. A degree of freedom on a corner of a domain
 * periodic in k directions is then mapped to the other 2^k - 1 corners. The partners
 * are stored in ascending order, so the degree of freedom with the smallest index
 * in a group is the only one with an index smaller than all of its partners.
 */
template<class Index>
void closePeriodicDofMap(std::unordered_map<Index, std::vector<Index>>& periodicDofMap)
{
    std::unordered_map<Index, std::vector<Index>> closedMap;
    for (const auto& entry : periodicDofMap)
    {
        if (closedMap.contains(entry.first))
            continue;

        std::vector<Index> group{ entry.first };
        for (std::size_t k = 0; k < group.size(); ++k)
            if (const auto it = periodicDofMap.find(group[k]); it != periodicDofMap.end())
                for (const auto periodicDofIdx : it->second)
                    if (std::ranges::find(group, periodicDofIdx) == group.end())
                        group.push_back(periodicDofIdx);

        std::ranges::sort(group);
        for (const auto dofIdx : group)
        {
            auto& periodicDofs = closedMap[dofIdx];
            for (const auto periodicDofIdx : group)
                if (periodicDofIdx != dofIdx)
                    periodicDofs.push_back(periodicDofIdx);
        }
    }

    periodicDofMap = std::move(closedMap);
}

} // end namespace Dumux::Detail

#endif
