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
public:
    using GlobalIndexSet = G;
    using LocalIndexSet = L;
    using GlobalIndexType = typename GlobalIndexSet::value_type;
    using LocalIndexType = typename LocalIndexSet::value_type;

    //! constructor fully defining the scv seed
    CCMpfaOLocalScvSeed(GlobalIndexSet&& globalScvfIndices,
                        LocalIndexSet&& localScvfIndices,
                        GlobalIndexType globalScvIndex)
    : globalScvIndex_(globalScvIndex),
      localScvfIndices_(std::move(localScvfIndices)),
      globalScvfIndices_(std::move(globalScvfIndices)) {}

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

/*!
 * \ingroup Mpfa
 * \brief Class for a sub-control volume face seed of the mpfa-o method.
 *
 * \param G the global index set type
 * \param L the local index set type
 */
template<typename G, typename L>
class CCMpfaOLocalScvfSeed
{
public:
    using GlobalIndexSet = G;
    using LocalIndexSet = L;
    using GlobalIndexType = typename GlobalIndexSet::value_type;
    using LocalIndexType = typename LocalIndexSet::value_type;

    //! constructor fully defining the scv face seed
    template<class SubControlVolumeFace>
    CCMpfaOLocalScvfSeed(const SubControlVolumeFace& scvf,
                         LocalIndexType insideLocalScvIndex,
                         LocalIndexSet&& outsideLocalScvIndices,
                         GlobalIndexSet&& outsideGlobalScvfIndices,
                         const MpfaFaceTypes faceType)
    : boundary_(scvf.boundary()),
      insideScvLocalIdx_(insideLocalScvIndex),
      outsideLocalScvIndices_(std::move(outsideLocalScvIndices)),
      insideScvGlobalIdx_(scvf.insideScvIdx()),
      outsideGlobalScvIndices_(scvf.outsideScvIndices()),
      insideScvfGlobalIdx_(scvf.index()),
      outsideGlobalScvfIndices_(std::move(outsideGlobalScvfIndices)),
      faceType_(faceType) {}

    //! Constructor when the outside indices are not known at this point
    template<class SubControlVolumeFace>
    CCMpfaOLocalScvfSeed(const SubControlVolumeFace& scvf,
                         LocalIndexType insideLocalScvIndex,
                         const MpfaFaceTypes faceType)
    : boundary_(scvf.boundary()),
      insideScvLocalIdx_(insideLocalScvIndex),
      insideScvGlobalIdx_(scvf.insideScvIdx()),
      outsideGlobalScvIndices_(scvf.outsideScvIndices()),
      insideScvfGlobalIdx_(scvf.index()),
      faceType_(faceType)
      {
        auto size = outsideGlobalScvIndices_.size();
        outsideLocalScvIndices_.reserve(size);
        outsideGlobalScvfIndices_.reserve(size);
      }

    const GlobalIndexType insideGlobalScvfIndex() const
    { return insideScvfGlobalIdx_; }

    const GlobalIndexSet& outsideGlobalScvfIndices() const
    { return outsideGlobalScvfIndices_; }

    const GlobalIndexType outsideGlobalScvfIndex() const
    {
        assert(outsideGlobalScvfIndices_.size() == 1 && "outside global scvf index not uniquely defined");
        return outsideGlobalScvfIndices_[0];
    }

    const LocalIndexType insideLocalScvIndex() const
    { return insideScvLocalIdx_; }

    const LocalIndexSet& outsideLocalScvIndices() const
    { return outsideLocalScvIndices_; }

    const GlobalIndexType insideGlobalScvIndex() const
    { return insideScvGlobalIdx_; }

    const GlobalIndexSet& outsideGlobalScvIndices() const
    { return outsideGlobalScvIndices_; }

    MpfaFaceTypes faceType() const
    { return faceType_; }

    bool boundary() const
    { return boundary_; }

    void addOutsideData(const GlobalIndexType outsideGlobalScvfIndex,
                        const LocalIndexType outsideLocalScvIndex)
    {
        outsideLocalScvIndices_.push_back(outsideLocalScvIndex);
        outsideGlobalScvfIndices_.push_back(outsideGlobalScvfIndex);
    }

    // for grids with dim < dimWorld, some outside indices might be doubled
    // we want to make the outside indices unique, but, the i-th outside global scvf face
    // should correspond to the j-th outside local scv.Therefore we apply the same operations on both containers
    void makeOutsideDataUnique()
    {
        for (auto scvIt = outsideLocalScvIndices_.begin(); scvIt != outsideLocalScvIndices_.end(); ++scvIt)
        {
            auto scvfIt = outsideGlobalScvfIndices_.begin();
            for (auto scvIt2 = scvIt+1; scvIt2 != outsideLocalScvIndices_.end(); ++scvIt2)
            {
                if (*scvIt2 == *scvIt)
                {
                    outsideLocalScvIndices_.erase(scvIt2);
                    outsideGlobalScvfIndices_.erase(scvfIt);
                    break;
                }
                ++scvfIt;
            }
        }
    }

private:
    bool boundary_;

    LocalIndexType insideScvLocalIdx_;
    LocalIndexSet outsideLocalScvIndices_;

    GlobalIndexType insideScvGlobalIdx_;
    GlobalIndexSet outsideGlobalScvIndices_;

    GlobalIndexType insideScvfGlobalIdx_;
    GlobalIndexSet outsideGlobalScvfIndices_;

    MpfaFaceTypes faceType_;
};
} // end namespace

#endif
