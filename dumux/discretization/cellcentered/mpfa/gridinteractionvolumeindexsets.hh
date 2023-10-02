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

#include "dualgridindexset.hh"

namespace Dumux {

template<class GridView>
class CCMpfaInteractionVolume
{
public:
    virtual ~CCMpfaInteractionVolume() = default;

    using GridIndex = typename GridView::IndexSet::IndexType;
    using GridIndexVisitor = std::function<void(const GridIndex&)>;

    void visitGridScvfIndices(const GridIndexVisitor& v) const
    { return visitGridScvfIndices_(v); }

private:
    virtual void visitGridScvfIndices_(const GridIndexVisitor&) const = 0;
};

template<typename GridView>
class CCMpfaInteractionVolumeFactory
{
public:
    virtual ~CCMpfaInteractionVolumeFactory() = default;

    using InteractionVolumePointer = std::unique_ptr<CCMpfaInteractionVolume<GridView>>;
    using DualGridNodalIndexSet = CCMpfaDualGridNodalIndexSet<GridView>;
    using InteractionVolumesVisitor = std::function<void(InteractionVolumePointer&&)>;

    void visitInteractionVolumesAt(const DualGridNodalIndexSet& ni,
                                   const InteractionVolumesVisitor& v) const
    { visitInteractionVolumesAt_(ni, v); }

private:
    virtual void visitInteractionVolumesAt_(const DualGridNodalIndexSet& ni,
                                            const InteractionVolumesVisitor&) const = 0;
};

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

    class DummyInteractionVolumeFactory : public CCMpfaInteractionVolumeFactory<GV>
    {
        class DummyInteractionVolume : public CCMpfaInteractionVolume<GV>
        {
        public:
            using typename CCMpfaInteractionVolume<GV>::GridIndexVisitor;
            DummyInteractionVolume(const CCMpfaDualGridNodalIndexSet<GV>& ni)
            : ni_{ni}
            {}

        private:
            void visitGridScvfIndices_(const GridIndexVisitor& v) const
            {
                for (const auto& scvfIndex : ni_.gridScvfIndices())
                    v(scvfIndex);
            }

            const CCMpfaDualGridNodalIndexSet<GV>& ni_;
        };

        using typename CCMpfaInteractionVolumeFactory<GV>::DualGridNodalIndexSet;
        using typename CCMpfaInteractionVolumeFactory<GV>::InteractionVolumesVisitor;
        void visitInteractionVolumesAt_(const DualGridNodalIndexSet& ni,
                                        const InteractionVolumesVisitor& v) const
        {
            v(std::make_unique<DummyInteractionVolume>(ni));
        }
    };

public:
    using GridGeometry = FVG;
    using PrimaryInteractionVolume = PI;
    using SecondaryInteractionVolume = SI;

    using DualGridIndexSet = CCMpfaDualGridIndexSet< GV >;
    using NodalIndexSet = typename DualGridIndexSet::NodalIndexSet;
    using GridIndexType = typename DualGridIndexSet::GridIndexType;

    void update(GridGeometry& gridGeometry,
                DualGridIndexSet&& dualGridIdSet)
    {
        DummyInteractionVolumeFactory f{};
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
                const CCMpfaInteractionVolumeFactory<GV>& factory)
    {
        dualGridIndexSet_ = std::make_unique<DualGridIndexSet>(std::move(dualGridIdSet));
        scvfIndexMap_.clear();
        scvfIndexMap_.resize(gridGeometry.numScvf());
        for (const auto& vertex : vertices(gridGeometry.gridView()))
            factory.visitInteractionVolumesAt((*dualGridIndexSet_)[gridGeometry.vertexMapper().index(vertex)], [&] (auto&& iv) {
                const auto ivIndex = interactionVolumes_.size();
                interactionVolumes_.emplace_back(std::move(iv));
                interactionVolumes_.back()->visitGridScvfIndices([&] (auto scvfIdx) {
                    scvfIndexMap_[scvfIdx] = ivIndex;
                });
            });
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
    std::vector<std::unique_ptr<CCMpfaInteractionVolume<GV>>> interactionVolumes_;
};

} // end namespace Dumux

#endif
