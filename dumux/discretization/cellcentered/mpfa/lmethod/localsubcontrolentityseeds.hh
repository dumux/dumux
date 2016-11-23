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
 * \brief Base class for sub control entity seeds of the mpfa-o method.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_L_LOCALSUBCONTROLENTITYSEEDS_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_L_LOCALSUBCONTROLENTITYSEEDS_HH

#include <dumux/implicit/cellcentered/mpfa/properties.hh>
#include <dumux/discretization/cellcentered/mpfa/facetypes.hh>

namespace Dumux
{
//! The scv seed class for scvs which could be the central scv in the actual interaction region
template<typename G, typename L>
class CCMpfaLCentralLocalScvSeed
{
    using GlobalIndexSet = G;
    using LocalIndexSet = L;
public:
    using GlobalIndexType = typename GlobalIndexSet::value_type;
    using LocalIndexType = typename LocalIndexSet::value_type;

    //! constructor fully defining the scv seed
    CCMpfaLCentralLocalScvSeed(GlobalIndexSet&& globalScvfIndices,
                               GlobalIndexType globalScvIndex,
                               LocalIndexType contiFaceLocalIdx)
    : contiFaceLocalIdx_(contiFaceLocalIdx),
      globalScvIndex_(globalScvIndex),
      globalScvfIndices_(std::move(globalScvfIndices)) {}

    LocalIndexType contiFaceLocalIdx() const
    { return contiFaceLocalIdx_; }

    GlobalIndexType globalIndex() const
    { return globalScvIndex_; }

    const GlobalIndexSet& globalScvfIndices() const
    { return globalScvfIndices_; }

private:
    LocalIndexType contiFaceLocalIdx_;
    GlobalIndexType globalScvIndex_;
    GlobalIndexSet globalScvfIndices_;
};

//! The scv seed class for scvs which could be the central scv in the actual interaction region
template<typename G>
class CCMpfaLOuterLocalScvSeed
{
    using GlobalIndexSet = G;
public:
    using GlobalIndexType = typename G::value_type;

    //! Constructor fully initializing the members
    CCMpfaLOuterLocalScvSeed(GlobalIndexType globalScvIndex,
                             GlobalIndexType globalScvfIndex)
    : scvIndexGlobal_(globalScvIndex),
      scvfIndexGlobal_(globalScvfIndex) {}

    //! Construct from central scv seed
    template<typename GI, typename LI>
    CCMpfaLOuterLocalScvSeed(const CCMpfaLCentralLocalScvSeed<GI, LI>& scvSeed)
    : scvIndexGlobal_(scvSeed.globalIndex()),
      scvfIndexGlobal_(scvSeed.globalScvfIndices()[scvSeed.contiFaceLocalIdx()]) {}


    const GlobalIndexType globalIndex() const
    { return scvIndexGlobal_; }

    const GlobalIndexType globalScvfIndex() const
    { return scvfIndexGlobal_; }

private:
    GlobalIndexType scvIndexGlobal_;
    GlobalIndexType scvfIndexGlobal_;
};

} // end namespace

#endif
