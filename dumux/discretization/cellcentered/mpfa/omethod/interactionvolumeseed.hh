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
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_O_INTERACTIONVOLUMESEED_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_O_INTERACTIONVOLUMESEED_HH

#include <dumux/implicit/cellcentered/mpfa/properties.hh>
#include "localsubcontrolentityseeds.hh"

namespace Dumux
{

/*!
 * \ingroup Mpfa
 * \brief Base class for the interaction volume seed of mpfa methods
 */
template<typename G, typename L>
class CCMpfaOInteractionVolumeSeed
{
    using GlobalIndexType = G;

public:
    using LocalScvSeed = CCMpfaOLocalScvSeed<G, L>;
    using LocalScvfSeed = CCMpfaOLocalScvfSeed<G, L>;

    CCMpfaOInteractionVolumeSeed(std::vector<LocalScvSeed>&& scvSeeds,
                                 std::vector<LocalScvfSeed>&& scvfSeeds,
                                 bool onBoundary)
    : onBoundary_(onBoundary),
      scvSeeds_(std::move(scvSeeds)),
      scvfSeeds_(std::move(scvfSeeds))
    {}

    bool onBoundary() const
    { return onBoundary_; }

    const std::vector<LocalScvSeed>& scvSeeds() const
    { return scvSeeds_; }

    const std::vector<LocalScvfSeed>& scvfSeeds() const
    { return scvfSeeds_; }

    std::vector<GlobalIndexType> globalScvIndices() const
    {
        std::vector<GlobalIndexType> globalIndices;
        globalIndices.reserve(scvSeeds().size());

        for (auto&& localScvSeed : scvSeeds())
            globalIndices.push_back(localScvSeed.globalIndex());

        return globalIndices;
    }

    std::vector<GlobalIndexType> globalScvfIndices() const
    {
        std::vector<GlobalIndexType> globalIndices;
        globalIndices.reserve(scvfSeeds().size() * 2);

        for (auto&& localScvfSeed : scvfSeeds())
            for (auto scvfIdxGlobal : localScvfSeed.globalScvfIndices())
                globalIndices.push_back(scvfIdxGlobal);

        globalIndices.shrink_to_fit();
        return globalIndices;
    }

private:
    bool onBoundary_;
    std::vector<LocalScvSeed> scvSeeds_;
    std::vector<LocalScvfSeed> scvfSeeds_;
};
} // end namespace

#endif
