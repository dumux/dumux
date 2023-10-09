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
#include <functional>

#include <dune/common/fmatrix.hh>

#include "dualgridindexset.hh"
#include "interactionvolume.hh"

// temporarily include and use o-method specific stuff
#include "mpfa_o.hh"

namespace Dumux {

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Class that holds all interaction volume index sets on a grid view.
 * \tparam FVG the finite volume grid geometry
 */
template<class FVG>
class CCMpfaGridInteractionVolumes
{
    using SubControlVolumeFace = typename FVG::SubControlVolumeFace;
    using GV = typename FVG::GridView;

public:
    using GridGeometry = FVG;

    using DualGridIndexSet = CCMpfaDualGridIndexSet<GV>;
    using NodalIndexSet = typename DualGridIndexSet::NodalIndexSet;
    using GridIndexType = typename DualGridIndexSet::GridIndexType;

    // for bw-compatibility
    void update(GridGeometry& gridGeometry,
                DualGridIndexSet&& dualGridIdSet)
    {
        CCMpfaO::InteractionVolumeFactory<FVG, double> f{};
        update(gridGeometry, std::move(dualGridIdSet), f);
    }

    /*!
     * \brief Construct all interaction volume index sets on the grid view
     *
     * \param gridGeometry The finite volume geometry on the grid view
     * \param dualGridIdSet The index sets of the dual grid on the grid view
     * \param ivFactory Factory class for constructing interaction volumes around vertices.
     */
    void update(GridGeometry& gridGeometry,
                DualGridIndexSet&& dualGridIdSet,
                const CCMpfaInteractionVolumeFactory<FVG, double>& factory)
    {
        dualGridIndexSet_ = std::make_unique<DualGridIndexSet>(std::move(dualGridIdSet));
        scvfIndexMap_.clear();
        scvfIndexMap_.resize(gridGeometry.numScvf());
        for (const auto& vertex : vertices(gridGeometry.gridView()))
            factory.visitInteractionVolumesAt((*dualGridIndexSet_)[gridGeometry.vertexMapper().index(vertex)], [&] (auto&& iv) {
                const auto ivIndex = interactionVolumes_.size();
                interactionVolumes_.emplace_back(std::move(iv));
                interactionVolumes_.back()->visitFluxGridScvfIndices([&] (auto scvfIdx) {
                    scvfIndexMap_[scvfIdx] = ivIndex;
                });
            });
    }

    //! Return the iv index set in which a given scvf is embedded in
    const auto& get(const SubControlVolumeFace& scvf) const
    { return at(scvf); }

    //! Return the iv index set in which a given scvf (index) is embedded in
    const auto& get(const GridIndexType scvfIdx) const
    { return at(scvfIdx); }

    //! Return the iv in which a given scvf is embedded in
    const CCMpfaInteractionVolume<FVG, double>& at(const SubControlVolumeFace& scvf) const
    { return at(scvf.index()); }

    //! Return the iv which a given scvf (index) is embedded in
    const CCMpfaInteractionVolume<FVG, double>& at(const GridIndexType scvfIdx) const
    { return *(interactionVolumes_[scvfIndexMap_[scvfIdx]]); }

    std::size_t size() const
    { return interactionVolumes_.size(); }

private:
    std::vector<GridIndexType> scvfIndexMap_;
    std::unique_ptr<DualGridIndexSet> dualGridIndexSet_;
    std::vector<std::unique_ptr<CCMpfaInteractionVolume<FVG, double>>> interactionVolumes_;
};

} // end namespace Dumux

#endif
