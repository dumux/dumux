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
 * \brief Class for the grid interaction volume index sets of mpfa schemes.
 */
#ifndef DUMUX_DISCRETIZATION_MPFA_O_GRIDINTERACTIONVOLUME_INDEXSETS_HH
#define DUMUX_DISCRETIZATION_MPFA_O_GRIDINTERACTIONVOLUME_INDEXSETS_HH

#include <dumux/discretization/cellcentered/mpfa/dualgridindexset.hh>

namespace Dumux
{
/*!
 * \ingroup Mpfa
 * \brief The grid interaction volume index sets class for the mpfa-o scheme.
 */
template<class TypeTag>
class CCMpfaGridInteractionVolumeIndexSets
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using PrimaryIV = typename GET_PROP_TYPE(TypeTag, PrimaryInteractionVolume);
    using PrimaryIVIndexSet = typename PrimaryIV::Traits::IndexSet;
    using SecondaryIV = typename GET_PROP_TYPE(TypeTag, SecondaryInteractionVolume);
    using SecondaryIVIndexSet = typename SecondaryIV::Traits::IndexSet;
    using Helper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    using DualGridIndexSet = CCMpfaDualGridIndexSet<TypeTag>;

    static const int dim = GridView::dimension;

public:

    void update(FVGridGeometry& fvGridGeometry, DualGridIndexSet&& dualGridIdSet)
    {
        dualGridIndexSet_ = std::make_unique<DualGridIndexSet>(std::move(dualGridIdSet));

        // clear containers
        primaryIVIndexSets_.clear();
        secondaryIVIndexSets_.clear();
        scvfIndexMap_.clear();

        // find out how many primary & secondary interaction volumes are needed
        numPrimaryIV_ = 0;
        numSecondaryIV_ = 0;
        for (const auto& vertex : vertices(fvGridGeometry.gridView()))
        {
            const auto vIdxGlobal = fvGridGeometry.vertexMapper().index(vertex);
            if (fvGridGeometry.vertexUsesPrimaryInteractionVolume(vIdxGlobal))
                numPrimaryIV_ += PrimaryIV::numInteractionVolumesAtVertex((*dualGridIndexSet_)[vIdxGlobal]);
            else
                numSecondaryIV_ += SecondaryIV::numInteractionVolumesAtVertex((*dualGridIndexSet_)[vIdxGlobal]);
        }

        // reserve memory
        primaryIVIndexSets_.reserve(numPrimaryIV_);
        secondaryIVIndexSets_.reserve(numSecondaryIV_);
        scvfIndexMap_.resize(fvGridGeometry.numScvf());

        // create interaction volume index sets around each vertex
        for (const auto& vertex : vertices(fvGridGeometry.gridView()))
        {
            const auto vIdxGlobal = fvGridGeometry.vertexMapper().index(vertex);
            if (fvGridGeometry.vertexUsesPrimaryInteractionVolume(vIdxGlobal))
                PrimaryIV::addInteractionVolumeIndexSets(primaryIVIndexSets_, scvfIndexMap_, (*dualGridIndexSet_)[vIdxGlobal]);
            else
                SecondaryIV::addInteractionVolumeIndexSets(secondaryIVIndexSets_, scvfIndexMap_, (*dualGridIndexSet_)[vIdxGlobal]);
        }
    }

    const PrimaryIVIndexSet& primaryIndexSet(const SubControlVolumeFace& scvf) const
    { return primaryIndexSet(scvf.index()); }

    const PrimaryIVIndexSet& primaryIndexSet(const IndexType scvfIdx) const
    { return primaryIVIndexSets_[scvfIndexMap_[scvfIdx]]; }

    const SecondaryIVIndexSet& secondaryIndexSet(const SubControlVolumeFace& scvf) const
    { return secondaryIndexSet(scvf.index()); }

    const SecondaryIVIndexSet& secondaryIndexSet(const IndexType scvfIdx) const
    { return secondaryIVIndexSets_[scvfIndexMap_[scvfIdx]]; }

    std::size_t numPrimaryInteractionVolumes() const
    { return numPrimaryIV_; }

    std::size_t numSecondaryInteractionVolumes() const
    { return numSecondaryIV_; }

private:
    std::vector<PrimaryIVIndexSet> primaryIVIndexSets_;
    std::vector<SecondaryIVIndexSet> secondaryIVIndexSets_;
    std::vector<IndexType> scvfIndexMap_;

    std::size_t numPrimaryIV_;
    std::size_t numSecondaryIV_;

    std::unique_ptr<DualGridIndexSet> dualGridIndexSet_;
};
} // end namespace


#endif
