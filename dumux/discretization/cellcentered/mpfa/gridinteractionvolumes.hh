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

namespace Dumux {

template<class GridGeometry, class Scalar>
class CCMpfaInteractionVolumeFactory
{
    using GridView = typename GridGeometry::GridView;

public:
    virtual ~CCMpfaInteractionVolumeFactory() = default;

    using DualGridNodalIndexSet = CCMpfaDualGridNodalIndexSet<GridView>;
    using InteractionVolume = CCMpfaInteractionVolume<GridGeometry, Scalar>;
    using InteractionVolumesVisitor = std::function<void(std::unique_ptr<InteractionVolume>&&)>;

    void visitInteractionVolumesAt(const DualGridNodalIndexSet& ni, const InteractionVolumesVisitor& v) const
    { visitInteractionVolumesAt_(ni, v); }

private:
    virtual void visitInteractionVolumesAt_(const DualGridNodalIndexSet&, const InteractionVolumesVisitor&) const = 0;
};

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

    class DummyTransmissibilities : public CCMpfaTransmissibilities<FVG, double>
    {
    public:
        using typename CCMpfaTransmissibilities<FVG, double>::SubControlVolumeFace;
        using typename CCMpfaTransmissibilities<FVG, double>::TensorAccessor;
        using typename CCMpfaTransmissibilities<FVG, double>::ForceAccessor;
        using typename CCMpfaTransmissibilities<FVG, double>::DofAccessor;

    private:
        int computeTransmissibilities_(const TensorAccessor&, const std::optional<ForceAccessor>&) override
        { DUNE_THROW(Dune::NotImplemented, "TransmissibilityFactory"); }

        double computeFlux_(const SubControlVolumeFace&, const DofAccessor&, int) const override
        { DUNE_THROW(Dune::NotImplemented, "TransmissibilityFactory"); }
    };

    class DummyInteractionVolumeFactory : public CCMpfaInteractionVolumeFactory<FVG, double>
    {
        class DummyInteractionVolume : public CCMpfaInteractionVolume<FVG, double>
        {
            using ParentType = CCMpfaInteractionVolume<FVG, double>;

        public:
            using typename CCMpfaInteractionVolume<FVG, double>::DirichletBoundaryPredicate;
            using typename CCMpfaInteractionVolume<FVG, double>::FVElementGeometry;
            using typename CCMpfaInteractionVolume<FVG, double>::GridIndexVisitor;

            DummyInteractionVolume(const CCMpfaDualGridNodalIndexSet<GV>& ni)
            : ni_{ni}
            {
                this->size_.emplace(typename ParentType::Size{
                    static_cast<unsigned int>(ni_.gridScvIndices().size()),  // numScvs
                    static_cast<unsigned int>(ni_.gridScvfIndices().size()), // numTotalScvfs
                    static_cast<unsigned int>(0),                            // numAuxiliaryScvfs
                    static_cast<unsigned int>(ni_.numBoundaryScvfs())        // numBoundaryScvfs
                });
            }

        private:
            bool isFluxScvf_(const SubControlVolumeFace&) const override
            { return true; }

            void visitGridScvIndices_(const GridIndexVisitor& v) const override
            {
                for (const auto& scvfIndex : ni_.gridScvIndices())
                    v(scvfIndex);
            }

            void visitGridScvfIndices_(const GridIndexVisitor& v) const override
            {
                for (const auto& scvfIndex : ni_.gridScvfIndices())
                    v(scvfIndex);
            }

            void visitFluxGridScvfIndices_(const GridIndexVisitor& v) const override
            { visitGridScvfIndices_(v); }

            std::unique_ptr<CCMpfaTransmissibilities<FVG, double>> transmissibilities_(
                const FVElementGeometry&,
                const DirichletBoundaryPredicate&
            ) const override
            { return std::make_unique<DummyTransmissibilities>(); }

            const CCMpfaDualGridNodalIndexSet<GV>& ni_;
        };

        using typename CCMpfaInteractionVolumeFactory<FVG, double>::DualGridNodalIndexSet;
        using typename CCMpfaInteractionVolumeFactory<FVG, double>::InteractionVolumesVisitor;
        void visitInteractionVolumesAt_(const DualGridNodalIndexSet& ni,
                                        const InteractionVolumesVisitor& v) const
        {
            v(std::make_unique<DummyInteractionVolume>(ni));
        }
    };

public:
    using GridGeometry = FVG;

    using DualGridIndexSet = CCMpfaDualGridIndexSet<GV>;
    using NodalIndexSet = typename DualGridIndexSet::NodalIndexSet;
    using GridIndexType = typename DualGridIndexSet::GridIndexType;

    // for bw-compatibility
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
    const NodalIndexSet& get(const SubControlVolumeFace& scvf) const
    { return get(scvf.index()); }

    //! Return the iv index set in which a given scvf (index) is embedded in
    const NodalIndexSet& get(const GridIndexType scvfIdx) const
    { return (*dualGridIndexSet_)[scvfIndexMap_[scvfIdx]]; }

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
