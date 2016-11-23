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
 * \brief Base class for interaction volume seeds of mpfa methods.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_L_INTERACTIONVOLUMESEED_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_L_INTERACTIONVOLUMESEED_HH

#include <dumux/implicit/cellcentered/mpfa/properties.hh>
#include "localsubcontrolentityseeds.hh"

namespace Dumux
{

/*!
 * \ingroup Mpfa
 * \brief Base class for the interaction volume seed of mpfa methods
 */
template<typename G, typename L>
class CCMpfaLInteractionVolumeSeed
{
    using GlobalIndexSet = G;
    using GlobalIndexType = typename G::value_type;
    using LocalIndexSet = L;
    using LocalIndexType = typename L::value_type;

public:
    using LocalScvSeed = CCMpfaLCentralLocalScvSeed<G, L>;
    using LocalOuterScvSeed = CCMpfaLOuterLocalScvSeed<G>;

    CCMpfaLInteractionVolumeSeed(std::vector<LocalScvSeed>&& scvSeeds,
                                 std::vector<LocalOuterScvSeed>&& outerScvSeeds,
                                 std::vector<GlobalIndexType>&& globalScvfIndices)
    : scvSeeds_(std::move(scvSeeds)),
      outerScvSeeds_(std::move(outerScvSeeds)),
      globalScvfIndices_(std::move(globalScvfIndices)) {}

    const std::vector<LocalScvSeed>& scvSeeds() const
    { return scvSeeds_; }

    const LocalScvSeed& scvSeed(const LocalIndexType idx) const
    { return scvSeeds_[idx]; }

    const std::vector<LocalOuterScvSeed>& outerScvSeeds() const
    { return outerScvSeeds_; }

    const LocalOuterScvSeed& outerScvSeed(const LocalIndexType idx) const
    { return outerScvSeeds_[idx]; }

    std::vector<GlobalIndexType> globalScvIndices() const
    {
        std::vector<GlobalIndexType> globalIndices;
        globalIndices.reserve(scvSeeds().size() + outerScvSeeds().size());

        for (auto&& localScvSeed : scvSeeds())
            globalIndices.push_back(localScvSeed.globalIndex());

        for (auto&& localScvSeed : outerScvSeeds())
            globalIndices.push_back(localScvSeed.globalIndex());

        // make the entries unique
        std::sort(globalIndices.begin(), globalIndices.end());
        globalIndices.erase(std::unique(globalIndices.begin(), globalIndices.end()), globalIndices.end());

        return globalIndices;
    }

    const std::vector<GlobalIndexType>& globalScvfIndices() const
    { return globalScvfIndices_; }

    bool isUnique() const
    { return scvSeeds_.size() == 1; }

private:
    std::vector<LocalScvSeed> scvSeeds_;
    std::vector<LocalOuterScvSeed> outerScvSeeds_;
    std::vector<GlobalIndexType> globalScvfIndices_;
};
} // end namespace

#endif
