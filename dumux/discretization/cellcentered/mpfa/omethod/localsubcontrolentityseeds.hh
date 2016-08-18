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
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_O_LOCALSUBCONTROLENTITYSEEDS_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_O_LOCALSUBCONTROLENTITYSEEDS_HH

#include <dumux/implicit/cellcentered/mpfa/properties.hh>
#include <dumux/discretization/cellcentered/mpfa/facetypes.hh>

namespace Dumux
{
template<typename G, typename L>
class CCMpfaOLocalScvSeed
{
    using GlobalIndexType = G;
    using LocalIndexType = L;
    using GlobalIndexSet = std::vector<GlobalIndexType>;
    using LocalIndexSet = std::vector<LocalIndexType>;

public:
    //! constructor fully defining the scv seed
    CCMpfaOLocalScvSeed(GlobalIndexSet&& globalScvfIndices,
                        LocalIndexSet&& localScvfIndices,
                        GlobalIndexType globalScvIndex)
    : globalScvIndex_(globalScvIndex),
      localScvfIndices_(std::move(localScvfIndices)),
      globalScvfIndices_(std::move(globalScvfIndices))
    {}

    //! Constructor when the local scvf indices are not known at this point
    CCMpfaOLocalScvSeed(GlobalIndexSet&& globalScvfIndices,
                        GlobalIndexType globalScvIndex)
    : globalScvIndex_(globalScvIndex),
      globalScvfIndices_(std::move(globalScvfIndices))
    {
        localScvfIndices_.resize(globalScvfIndices_.size(), -1);
    }

    GlobalIndexType globalIndex() const
    { return globalScvIndex_; }

    const LocalIndexSet& localScvfIndices() const
    { return localScvfIndices_; }

    const GlobalIndexSet& globalScvfIndices() const
    { return globalScvfIndices_; }

    void setLocalScvfIndex(int coordDir, LocalIndexType localScvfIdx)
    {
        localScvfIndices_[coordDir] = localScvfIdx;
    }

private:
    GlobalIndexType globalScvIndex_;
    LocalIndexSet localScvfIndices_;
    GlobalIndexSet globalScvfIndices_;
};


template<typename G, typename L>
class CCMpfaOLocalScvfSeed
{
    using GlobalIndexType = G;
    using LocalIndexType = L;
    using GlobalIndexSet = std::vector<GlobalIndexType>;
    using LocalIndexSet = std::vector<LocalIndexType>;

public:
    template<class SubControlVolumeFace>
    CCMpfaOLocalScvfSeed(const SubControlVolumeFace& scvf,
                         LocalIndexSet&& scvIndicesLocal,
                         GlobalIndexSet&& scvfIndicesGlobal,
                         const MpfaFaceTypes faceType)
    : boundary_(scvf.boundary()),
      scvIndicesGlobal_({scvf.insideScvIdx(), scvf.outsideScvIdx()}),
      scvfIndicesGlobal_(std::move(scvfIndicesGlobal)),
      scvIndicesLocal_(std::move(scvIndicesLocal)),
      faceType_(faceType)
    {}

    const GlobalIndexSet& globalScvfIndices() const
    { return scvfIndicesGlobal_; }

    const LocalIndexSet& localScvIndices() const
    { return scvIndicesLocal_; }

    const GlobalIndexSet& globalScvIndices() const
    { return scvIndicesGlobal_; }

    LocalIndexType insideLocalScvIndex() const
    { return scvIndicesLocal_[0]; }

    GlobalIndexType insideGlobalScvfIndex() const
    { return scvfIndicesGlobal_[0]; }

    MpfaFaceTypes faceType() const
    { return faceType_; }

    bool boundary() const
    { return boundary_; }

private:
    bool boundary_;
    GlobalIndexSet scvIndicesGlobal_;
    GlobalIndexSet scvfIndicesGlobal_;
    LocalIndexSet scvIndicesLocal_;
    MpfaFaceTypes faceType_;
};
} // end namespace

#endif
