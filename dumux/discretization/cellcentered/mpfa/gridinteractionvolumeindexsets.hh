// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup CCMpfaDiscretization
 * \brief Class for the grid interaction volume index sets of mpfa schemes.
 */
#ifndef DUMUX_DISCRETIZATION_MPFA_O_GRIDINTERACTIONVOLUME_INDEXSETS_HH
#define DUMUX_DISCRETIZATION_MPFA_O_GRIDINTERACTIONVOLUME_INDEXSETS_HH

#include <memory>

#include "dualgridindexset.hh"

namespace Dumux {

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Class that holds all interaction volume index sets on a grid view.
 *
 * \tparam FVG the finite volume grid geometry
 * \tparam NI the type used for nodal index sets
 * \tparam PI primary interaction volume type
 * \tparam SI secondary interaction volume type
 */
template<class FVG, class NI, class PI, class SI = PI>
class CCMpfaGridInteractionVolumeIndexSets
{
    using SubControlVolumeFace = typename FVG::SubControlVolumeFace;
    using PrimaryIVIndexSet = typename PI::Traits::IndexSet;
    using SecondaryIVIndexSet = typename SI::Traits::IndexSet;

public:
    using GridGeometry = FVG;
    using PrimaryInteractionVolume = PI;
    using SecondaryInteractionVolume = SI;

    using GridIndexType = typename NI::GridIndexType;
    using DualGridIndexSet = CCMpfaDualGridIndexSet< NI >;

    /*!
     * \brief Construct all interaction volume index sets on the grid view
     *
     * \param gridGeometry The finite volume geometry on the grid view
     * \param dualGridIdSet The index sets of the dual grid on the grid view
     */
    void update(GridGeometry& gridGeometry, DualGridIndexSet&& dualGridIdSet)
    {
        dualGridIndexSet_ = std::make_unique<DualGridIndexSet>(std::move(dualGridIdSet));

        // clear containers
        primaryIVIndexSets_.clear();
        secondaryIVIndexSets_.clear();
        scvfIndexMap_.clear();

        // find out how many primary & secondary interaction volumes are needed
        numPrimaryIV_ = 0;
        numSecondaryIV_ = 0;
        for (const auto& vertex : vertices(gridGeometry.gridView()))
        {
            const auto vIdxGlobal = gridGeometry.vertexMapper().index(vertex);
            if (!gridGeometry.vertexUsesSecondaryInteractionVolume(vIdxGlobal))
                numPrimaryIV_ += PrimaryInteractionVolume::numIVAtVertex((*dualGridIndexSet_)[vIdxGlobal]);
            else
                numSecondaryIV_ += SecondaryInteractionVolume::numIVAtVertex((*dualGridIndexSet_)[vIdxGlobal]);
        }

        // reserve memory
        primaryIVIndexSets_.reserve(numPrimaryIV_);
        secondaryIVIndexSets_.reserve(numSecondaryIV_);
        scvfIndexMap_.resize(gridGeometry.numScvf());

        // create interaction volume index sets around each vertex
        for (const auto& vertex : vertices(gridGeometry.gridView()))
        {
            const auto vIdxGlobal = gridGeometry.vertexMapper().index(vertex);
            if (!gridGeometry.vertexUsesSecondaryInteractionVolume(vIdxGlobal))
                PrimaryInteractionVolume::addIVIndexSets(primaryIVIndexSets_,
                                                         scvfIndexMap_,
                                                         (*dualGridIndexSet_)[vIdxGlobal],
                                                         gridGeometry.flipScvfIndexSet());
            else
                SecondaryInteractionVolume::addIVIndexSets(secondaryIVIndexSets_,
                                                           scvfIndexMap_,
                                                           (*dualGridIndexSet_)[vIdxGlobal],
                                                           gridGeometry.flipScvfIndexSet());
        }
    }

    //! Return the iv index set in which a given scvf is embedded in
    const PrimaryIVIndexSet& primaryIndexSet(const SubControlVolumeFace& scvf) const
    { return primaryIndexSet(scvf.index()); }

    //! Return the iv index set in which a given scvf (index) is embedded in
    const PrimaryIVIndexSet& primaryIndexSet(const GridIndexType scvfIdx) const
    { return primaryIVIndexSets_[scvfIndexMap_[scvfIdx]]; }

    //! Return the iv index set in which a given scvf is embedded in
    const SecondaryIVIndexSet& secondaryIndexSet(const SubControlVolumeFace& scvf) const
    { return secondaryIndexSet(scvf.index()); }

    //! Return the iv index set in which a given scvf (index) is embedded in
    const SecondaryIVIndexSet& secondaryIndexSet(const GridIndexType scvfIdx) const
    { return secondaryIVIndexSets_[scvfIndexMap_[scvfIdx]]; }

    //! Returns number of primary/secondary interaction volumes on the grid view
    std::size_t numPrimaryInteractionVolumes() const { return numPrimaryIV_; }
    std::size_t numSecondaryInteractionVolumes() const { return numSecondaryIV_; }

private:
    std::vector<PrimaryIVIndexSet> primaryIVIndexSets_;
    std::vector<SecondaryIVIndexSet> secondaryIVIndexSets_;
    std::vector<GridIndexType> scvfIndexMap_;

    std::size_t numPrimaryIV_;
    std::size_t numSecondaryIV_;

    std::unique_ptr<DualGridIndexSet> dualGridIndexSet_;
};

} // end namespace Dumux

#endif
