// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MultiDomain
 * \ingroup MortarCoupling
 * \brief Base problem for subdomain problems in mortar-coupling models.
 */
#ifndef DUMUX_MORTAR_SUBDOMAIN_PROBLEM_BASE_HH
#define DUMUX_MORTAR_SUBDOMAIN_PROBLEM_BASE_HH

#include <ranges>
#include <vector>
#include <optional>
#include <unordered_map>

#include <dune/grid/common/mcmgmapper.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/io/grid/facetgridmanager.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/traceoperator.hh>
#include <dumux/discretization/facetgridmapper.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux::Mortar {

#ifndef DOXYGEN
namespace Detail {

template<typename GV>
class TraceGridGeometryFacade
{
 public:
    using GridView = GV;

    TraceGridGeometryFacade(const GridView& gv)
    : gridView_{gv}
    , elementMapper_{gv, Dune::mcmgElementLayout()}
    , vertexMapper_{gv, Dune::mcmgVertexLayout()}
    {}

    const GridView& gridView() const { return gridView_; }
    const auto& elementMapper() const { return elementMapper_; }
    const auto& vertexMapper() const { return vertexMapper_; }

 private:
    GridView gridView_;
    Dune::MultipleCodimMultipleGeomTypeMapper<GridView> elementMapper_;
    Dune::MultipleCodimMultipleGeomTypeMapper<GridView> vertexMapper_;
};

}  // namespace Detail
#endif  // DOXYGEN

/*!
 * \ingroup MultiDomain
 * \ingroup MortarCoupling
 * \brief Base problem for finite-volume subdomain problems in mortar-coupling models.
 *        Provides default implementations for some of the interfaces
 *        called by the default subdomain solver.
 */
template<typename GridGeometry,
         typename MortarGrid,
         typename MortarSolution>
class SubDomainFVProblemBase
{
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;

    using TraceGridManager = FacetGridManager<typename GridView::Grid, MortarGrid>;
    using TraceGridMapper = FVFacetGridMapper<typename MortarGrid::LeafGridView, GridGeometry>;
    using TraceGridGeometry = Detail::TraceGridGeometryFacade<typename MortarGrid::LeafGridView>;
    using TraceOperator = FVTraceOperator<TraceGridGeometry, GridGeometry>;

    static constexpr bool isCVFE = DiscretizationMethods::isCVFE<typename GridGeometry::DiscretizationMethod>;

    struct EntityCouplingMap {
        std::size_t sceIndex;
        std::size_t mortarId;
        std::size_t traceDofIndex;
    };

public:
    using TraceSolutionVector = MortarSolution;

    SubDomainFVProblemBase(std::shared_ptr<const GridGeometry> gridGeometry)
    : gridGeometry_{std::move(gridGeometry)}
    {}

    //! Register a trace with the mortar domain that has the given id
    void registerMortarTrace(std::shared_ptr<const TraceGridManager> traceGridManager, std::size_t id)
    {
        const auto& traceGridView = traceGridManager->grid().leafGridView();
        auto traceGridMapper = std::make_shared<TraceGridMapper>(traceGridView, gridGeometry_);
        auto traceGridGeometry = std::make_shared<TraceGridGeometry>(traceGridView);

        std::vector<std::size_t> domainToTraceVertexIndexMap;
        if constexpr (isCVFE)
        {
            domainToTraceVertexIndexMap.resize(gridGeometry_->vertexMapper().size());
            for (const auto& v : vertices(traceGridView))
                domainToTraceVertexIndexMap[gridGeometry_->vertexMapper().index(traceGridManager->hostGridVertex(v))]
                    = traceGridGeometry->vertexMapper().index(v);
        }

        coupledEntities_.resize(gridGeometry_->gridView().size(0));
        for (const auto& mortarElement : elements(traceGridView))
            for (const auto& element : traceGridMapper->domainElementsAdjacentTo(mortarElement))
            {
                const auto eIdx = gridGeometry_->elementMapper().index(element);
                const auto fvGeometry = localView(*gridGeometry_).bindElement(element);
                if constexpr(isCVFE)
                    for (const auto scvfIdx : traceGridMapper->domainScvfsAdjacentTo(mortarElement, element))
                    {
                        const auto& scv = fvGeometry.scv(fvGeometry.scvf(scvfIdx).insideScvIdx());
                        coupledEntities_[eIdx].push_back({
                            .sceIndex = subControlEntityIndex_(scv),
                            .mortarId = id,
                            .traceDofIndex = domainToTraceVertexIndexMap.at(scv.dofIndex())
                        });
                    }
                else
                    for (const auto scvfIdx : traceGridMapper->domainScvfsAdjacentTo(mortarElement, element))
                        coupledEntities_[eIdx].push_back({
                            .sceIndex = subControlEntityIndex_(fvGeometry.scvf(scvfIdx)),
                            .mortarId = id,
                            .traceDofIndex = traceGridGeometry->elementMapper().index(mortarElement)
                        });
            }
        mortarIdToTraceOperator_.emplace(std::make_pair(
            id,
            TraceOperator{traceGridGeometry, traceGridMapper, [&] (const auto& traceVertex) {
                return traceGridManager->hostGridVertex(traceVertex);
            }}
        ));
    }

    //! Return true if the given element/sub-control volume face coincides with
    bool isOnMortarBoundary(const Element& element, const SubControlVolumeFace& scvf) const requires(!isCVFE)
    { return getCouplingMap_(element, scvf).has_value(); }

    //! Return true if the given element/sub-control volume coincides with
    bool isOnMortarBoundary(const Element& element, const SubControlVolume& scv) const requires(isCVFE)
    { return getCouplingMap_(element, scv).has_value(); }

    //! Set the mortar variables to be used on the trace with the mortar domain with the given id
    void setTraceVariables(std::size_t mortarId, TraceSolutionVector trace)
    { mortarBoundaryConditions_[mortarId] = std::move(trace); }

    //! Return the trace operator for the part of the trace overlapping with the mortar domain with the given id
    const auto& getTraceOperator(const std::size_t mortarId) const
    { return mortarIdToTraceOperator_.at(mortarId); }

    //! Return the trace variables associated with the given boundary scvf
    const auto& getTraceVariables(const Element& element, const SubControlVolumeFace& scvf) const requires(!isCVFE)
    {
        const auto& map = getCouplingMap_(element, scvf);
        if (!map.has_value())
            DUNE_THROW(Dune::InvalidStateException, "Given scvf does not overlap with a mortar domain.");
        return mortarBoundaryConditions_.at(map->mortarId)[map->traceDofIndex];
    }

    //! Return the trace variables associated with the given boundary scv
    auto getTraceVariables(const Element& element, const SubControlVolume& scv) const requires(isCVFE)
    {
        typename TraceSolutionVector::value_type result(0);
        int count = 0;
        visitCouplingMaps_(element, scv, [&] (const auto& entry) {
            result += mortarBoundaryConditions_.at(entry.mortarId)[entry.traceDofIndex];
            count++;
        });
        if (count == 0)
            DUNE_THROW(Dune::InvalidStateException, "Given scv does not overlap with a mortar domain.");
        result /= count;
        return result;
    }

private:
    template<typename SubControlEntity>
    std::optional<EntityCouplingMap> getCouplingMap_(const Element& element, const SubControlEntity& sce) const
    {
        if (coupledEntities_.empty())
            return {};

        const auto eIdx = gridGeometry_->elementMapper().index(element);
        const auto& map = coupledEntities_.at(eIdx);
        auto it = std::ranges::find_if(map, [&] (const auto& entry) { return subControlEntityIndex_(sce) == entry.sceIndex; });
        if (it == std::ranges::end(map))
            return {};
        return *it;
    }

    template<typename Visitor>
    void visitCouplingMaps_(const Element& element, const SubControlVolume& scv, Visitor&& v) const
    {
        const auto eIdx = gridGeometry_->elementMapper().index(element);
        std::ranges::for_each(coupledEntities_.at(eIdx) | std::views::filter([&] (const auto& entry) {
            return subControlEntityIndex_(scv) == entry.sceIndex;
        }), v);
    }

    std::size_t subControlEntityIndex_(const SubControlVolumeFace& scvf) const { return scvf.index(); }
    std::size_t subControlEntityIndex_(const SubControlVolume& scv) const { return scv.dofIndex(); }

    std::shared_ptr<const GridGeometry> gridGeometry_;
    std::vector<std::vector<EntityCouplingMap>> coupledEntities_;
    std::unordered_map<std::size_t, TraceOperator> mortarIdToTraceOperator_;
    std::unordered_map<std::size_t, TraceSolutionVector> mortarBoundaryConditions_;
};

} // end namespace Dumux::Mortar

#endif
