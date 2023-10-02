// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
 * \tparam PI primary interaction volume type
 * \tparam SI secondary interaction volume type
 */
template<class FVG, class PI, class SI = PI>
class CCMpfaGridInteractionVolumeIndexSets
{
    using SubControlVolumeFace = typename FVG::SubControlVolumeFace;
    using GV = typename FVG::GridView;

public:
    using GridGeometry = FVG;
    using PrimaryInteractionVolume = PI;
    using SecondaryInteractionVolume = SI;

    using DualGridIndexSet = CCMpfaDualGridIndexSet< GV >;
    using NodalIndexSet = typename DualGridIndexSet::NodalIndexSet;
    using GridIndexType = typename DualGridIndexSet::GridIndexType;

    /*!
     * \brief Construct all interaction volume index sets on the grid view
     *
     * \param gridGeometry The finite volume geometry on the grid view
     * \param dualGridIdSet The index sets of the dual grid on the grid view
     */
    void update(GridGeometry& gridGeometry, DualGridIndexSet&& dualGridIdSet)
    {
        dualGridIndexSet_ = std::make_unique<DualGridIndexSet>(std::move(dualGridIdSet));
        scvfIndexMap_.clear();
        scvfIndexMap_.resize(gridGeometry.numScvf());
        for (const auto& vertex : vertices(gridGeometry.gridView()))
        {
            const auto vIdxGlobal = gridGeometry.vertexMapper().index(vertex);
            for (const auto scvfIdx : (*dualGridIndexSet_)[vIdxGlobal].gridScvfIndices())
                scvfIndexMap_[scvfIdx] = vIdxGlobal;
        }
    }

    //! Return the iv index set in which a given scvf is embedded in
    const NodalIndexSet& get(const SubControlVolumeFace& scvf) const
    { return get(scvf.index()); }

    //! Return the iv index set in which a given scvf (index) is embedded in
    const NodalIndexSet& get(const GridIndexType scvfIdx) const
    { return (*dualGridIndexSet_)[scvfIndexMap_[scvfIdx]]; }

    std::size_t numPrimaryInteractionVolumes() const { return dualGridIndexSet_->size(); }
    std::size_t numSecondaryInteractionVolumes() const { return 0; }

private:
    std::vector<GridIndexType> scvfIndexMap_;
    std::unique_ptr<DualGridIndexSet> dualGridIndexSet_;
};

} // end namespace Dumux

#endif
