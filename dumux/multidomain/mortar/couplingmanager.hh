// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MortarCoupling
 * \brief Coupling manager of a subdomain in a mortar-coupled model.
 *
 * SubDomainCouplingManager is the name subdomain problems compose. It resolves to the
 * implementation for the trace kind of the subdomain, FacetTraceCouplingManager for bulk
 * subdomains, whose trace is a grid of boundary facets. Each implementation comes with
 * overloads of the free evaluation functions isOnMortarBoundary, traceVariablesAt and
 * traceValues, so code generic in the manager calls those.
 */
#ifndef DUMUX_MULTIDOMAIN_MORTAR_COUPLING_MANAGER_HH
#define DUMUX_MULTIDOMAIN_MORTAR_COUPLING_MANAGER_HH

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/reservedvector.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <dune/localfunctions/lagrange/lagrangelfecache.hh>

#include <dumux/common/concepts/field_.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/facetgridmapper.hh>
#include <dumux/geometry/intersectingentities.hh>

#include "couplingmode.hh"
#include "trace.hh"

namespace Dumux::Mortar {

/*!
 * \ingroup MortarCoupling
 * \brief The sub-entities of a subdomain a trace is assembled from: the boundary
 *        sub-control volume faces or the boundary faces of a bulk subdomain, or the
 *        sub-control volumes of the pores of a network.
 */
enum class TraceEntity { subControlVolumeFace, boundaryFace, subControlVolume };

/*!
 * \ingroup MortarCoupling
 * \brief The coupling of a bulk subdomain to its mortar domains through a facet trace: which
 *        sub-control entities lie on which trace cell, the mortar data imposed on the traces,
 *        the mode in which that data enters the subdomain problem, and the evaluation of the
 *        coupling conditions.
 *
 * There is one such manager per subdomain, built by the subdomain solver and composed by the
 * subdomain problem, which asks it whether an entity is coupled and what the mortar data
 * there is. Its evaluation members forward to the free functions of this file, which are
 * the implementation and serve code that is generic in the manager. Problems name the
 * manager through SubDomainCouplingManager, which selects the implementation matching the
 * trace kind of the subdomain.
 *
 * The data imposed on a trace is given either per trace cell or per trace vertex, as the
 * order of the trace data declares; in the latter case the trace function is linear within
 * each trace cell. Reads at faces and positions evaluate that function, so a problem reads
 * the same thing whatever the layout.
 */
template<typename GG, typename MortarGrid, typename MortarSolution>
class FacetTraceCouplingManager
{
    using GridView = typename GG::GridView;
    using TraceGridMapper = FacetGridMapper<typename MortarGrid::LeafGridView, GG>;

public:
    static constexpr bool isCVFE = DiscretizationMethods::isCVFE<typename GG::DiscretizationMethod>;

    using GridGeometry = GG;
    using Trace = Mortar::Trace<GG, MortarGrid>;
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using TraceSolutionVector = MortarSolution;
    using TraceValue = typename MortarSolution::value_type;
    using Scalar = typename MortarSolution::field_type;
    using TraceCellGeometry = Dune::MultiLinearGeometry<typename GlobalPosition::value_type, MortarGrid::dimension, GlobalPosition::dimension>;

    //! couples a sub-control entity (a face for cell-centred, a degree of freedom for
    //! control-volume finite element schemes) to a trace degree of freedom
    struct EntityCouplingMap {
        std::size_t sceIndex;
        std::size_t mortarId;
        std::size_t traceDofIndex;
        bool operator==(const EntityCouplingMap&) const = default;
    };

    //! couples a sub-control volume face to the trace cell it lies in
    struct ScvfCouplingMap {
        std::size_t scvfIndex;
        std::size_t mortarId;
        std::size_t traceCellIndex;
    };

    //! couples a boundary face of an element, by the index of its intersection, to the
    //! trace cell it lies in
    struct BoundaryFaceCouplingMap {
        std::size_t intersectionIndex;
        std::size_t mortarId;
        std::size_t traceCellIndex;
    };

    //! The manager of the subdomain with the given grid geometry, with no trace registered yet
    explicit FacetTraceCouplingManager(std::shared_ptr<const GridGeometry> gridGeometry)
    : gridGeometry_{std::move(gridGeometry)}
    {}

    //! Register the trace shared with the mortar domain that has the given id
    void registerTrace(std::shared_ptr<const Trace> trace, std::size_t id)
    {
        const auto& traceGridView = trace->gridView();
        const auto& traceCells = *trace->entitySet();
        const TraceGridMapper traceGridMapper(*trace, gridGeometry_);

        // keyed by the discretization's dof indices, which for hybrid schemes interleave
        // vertex and non-vertex degrees of freedom; entries that are no trace vertex stay
        // at the sentinel
        constexpr auto noTraceVertex = std::numeric_limits<std::size_t>::max();
        std::vector<std::size_t> domainDofToTraceVertexIndexMap;
        if constexpr (isCVFE)
        {
            const auto domainDofIndexOf = [&] (const auto& vertex) -> std::size_t {
                if constexpr (requires { gridGeometry_->dofMapper().index(vertex); })
                    return gridGeometry_->dofMapper().index(vertex);
                else
                    return gridGeometry_->vertexMapper().index(vertex);
            };
            domainDofToTraceVertexIndexMap.assign(gridGeometry_->numDofs(), noTraceVertex);
            for (const auto& v : vertices(traceGridView))
                domainDofToTraceVertexIndexMap[domainDofIndexOf(trace->hostGridVertex(v))]
                    = traceGridView.indexSet().index(v);
        }

        coupledEntities_.resize(gridGeometry_->gridView().size(0));
        coupledScvfs_.resize(gridGeometry_->gridView().size(0));
        coupledBoundaryFaces_.resize(gridGeometry_->gridView().size(0));

        // the cell geometries and vertices, in the order of the entity set's indices
        std::vector<std::vector<GlobalPosition>> cellCorners(traceCells.size());
        std::vector<Dune::GeometryType> cellTypes(traceCells.size());
        auto& cellVertices = traceCellVertices_[id];
        cellVertices.clear();
        cellVertices.resize(traceCells.size());
        for (const auto& traceCell : traceCells)
        {
            const auto traceCellIndex = traceCells.index(traceCell);
            const auto traceCellGeometry = traceCell.geometry();
            cellTypes[traceCellIndex] = traceCellGeometry.type();
            cellCorners[traceCellIndex].resize(traceCellGeometry.corners());
            for (int c = 0; c < traceCellGeometry.corners(); ++c)
            {
                cellCorners[traceCellIndex][c] = traceCellGeometry.corner(c);
                cellVertices[traceCellIndex].push_back(traceGridView.indexSet().subIndex(traceCell, c, MortarGrid::dimension));
            }
            // the cache creates a finite element on first use, which must not happen during
            // a concurrent assembly
            feCache_.get(traceCellGeometry.type());

            for (const auto& element : traceGridMapper.domainElementsAdjacentTo(traceCell))
            {
                const auto eIdx = gridGeometry_->elementMapper().index(element);
                const auto fvGeometry = localView(*gridGeometry_).bindElement(element);
                const auto face = traceGridMapper.boundaryFace(fvGeometry, traceCell);
                coupledBoundaryFaces_[eIdx].push_back({
                    .intersectionIndex = static_cast<std::size_t>(face.intersectionIndex()),
                    .mortarId = id,
                    .traceCellIndex = traceCellIndex
                });
                for (const auto& scvf : scvfs(fvGeometry, face))
                    coupledScvfs_[eIdx].push_back({
                        .scvfIndex = scvf.index(),
                        .mortarId = id,
                        .traceCellIndex = traceCellIndex
                    });

                if constexpr (isCVFE)
                {
                    // degrees of freedom without a trace vertex take part only in the
                    // face-wise coupling maps; one shared by two trace cells is recorded once
                    for (const auto& localDof : localDofs(fvGeometry, face))
                    {
                        const auto traceDofIndex = domainDofToTraceVertexIndexMap[localDof.dofIndex()];
                        if (traceDofIndex == noTraceVertex)
                            continue;
                        const EntityCouplingMap entry{static_cast<std::size_t>(localDof.dofIndex()), id, traceDofIndex};
                        if (std::ranges::find(coupledEntities_[eIdx], entry) == coupledEntities_[eIdx].end())
                            coupledEntities_[eIdx].push_back(entry);
                    }
                }
                else
                    for (const auto& scvf : scvfs(fvGeometry, face))
                        coupledEntities_[eIdx].push_back({
                            .sceIndex = static_cast<std::size_t>(scvf.index()),
                            .mortarId = id,
                            .traceDofIndex = traceCellIndex
                        });
            }
        }
        auto& cellGeometries = traceCellGeometries_[id];
        cellGeometries.clear();
        cellGeometries.reserve(traceCells.size());
        for (std::size_t i = 0; i < traceCells.size(); ++i)
            cellGeometries.emplace_back(cellTypes[i], std::move(cellCorners[i]));
        traces_[id] = std::move(trace);
        isFloating_.reset();
    }

    /*!
     * \brief Set the mortar data to be imposed on the trace shared with the given mortar
     *        domain: one entry per trace cell or per trace vertex, as the order of the trace
     *        data declares.
     */
    void setTraceVariables(std::size_t mortarId, TraceSolutionVector trace)
    {
        if (trace.size() != numTraceDofs(mortarId))
            DUNE_THROW(Dune::InvalidStateException, "Trace data for mortar " << mortarId << " has " << trace.size()
                       << " entries, but the trace takes " << numTraceDofs(mortarId)
                       << " entries for data of order " << traceDataOrder_);
        traceData_[mortarId] = std::move(trace);
        ++dataVersion_;
    }

    /*!
     * \brief The order of the trace data this manager reads: zero for one value per trace
     *        cell, one for one value per trace vertex, the trace function being linear
     *        within each trace cell.
     */
    std::size_t traceDataOrder() const
    { return traceDataOrder_; }

    //! Set the order of the trace data, before any data is imposed
    void setTraceDataOrder(std::size_t order)
    {
        if (order > 1)
            DUNE_THROW(Dune::NotImplemented, "Trace data of order " << order);
        traceDataOrder_ = order;
    }

    //! The number of entries the data on the trace shared with the given mortar domain has
    std::size_t numTraceDofs(std::size_t mortarId) const
    {
        const auto it = traces_.find(mortarId);
        if (it == traces_.end())
            DUNE_THROW(Dune::InvalidStateException, "No trace is registered for mortar " << mortarId);
        return traceDataOrder_ == 0 ? it->second->size() : it->second->gridView().size(MortarGrid::dimension);
    }

    //! A counter that changes whenever mortar data is set, for consumers caching data
    //! derived from it
    std::size_t dataVersion() const
    { return dataVersion_; }

    //! Set the mode in which the mortar data enters the subdomain problem
    void setCouplingMode(CouplingMode mode)
    {
        couplingMode_ = mode;
        // cache the flag here, in a single-threaded context, so that concurrent
        // assembly threads only ever read it
        if (mode == CouplingMode::natural)
            isFloating();
    }

    //! Return the mode in which the mortar data enters the subdomain problem
    CouplingMode couplingMode() const
    { return couplingMode_; }

    //! Return true if the entire boundary of the subdomain couples to mortar domains
    bool isFloating() const
    {
        if (!isFloating_.has_value())
        {
            isFloating_ = true;
            auto fvGeometry = localView(*gridGeometry_);
            for (const auto& element : elements(gridGeometry_->gridView()))
            {
                if (!element.hasBoundaryIntersections())
                    continue;
                fvGeometry.bindElement(element);
                const auto eIdx = gridGeometry_->elementMapper().index(element);
                for (const auto& scvf : scvfs(fvGeometry))
                {
                    if (!scvf.boundary())
                        continue;
                    if (std::ranges::none_of(faceCouplingsOf(eIdx),
                                             [&] (const auto& e) { return e.scvfIndex == scvf.index(); }))
                    {
                        isFloating_ = false;
                        return *isFloating_;
                    }
                }
            }
        }
        return *isFloating_;
    }

    //! Return the grid geometry of the subdomain
    const GridGeometry& gridGeometry() const
    { return *gridGeometry_; }

    //! Return the number of cells of the trace shared with the given mortar domain
    std::size_t numTraceCells(std::size_t mortarId) const
    { return traces_.at(mortarId)->size(); }

    //! Return the geometry of the given cell of the trace shared with the given mortar domain
    const TraceCellGeometry& traceCellGeometry(std::size_t mortarId, std::size_t traceCellIndex) const
    { return traceCellGeometries_.at(mortarId)[traceCellIndex]; }

    //! Return true if the given cell of the trace shared with the given mortar domain contains the position
    bool traceCellContains(std::size_t mortarId, std::size_t traceCellIndex, const GlobalPosition& globalPos) const
    {
        using ctype = typename GlobalPosition::value_type;
        static constexpr int traceDim = MortarGrid::dimension;
        const auto& geometry = traceCellGeometry(mortarId, traceCellIndex);
        const auto local = geometry.local(globalPos);
        const auto tolerance = 1e-7*std::pow(geometry.volume(), 1.0/traceDim);
        return Dune::referenceElement<ctype, traceDim>(geometry.type()).checkInside(local)
            && (geometry.global(local) - globalPos).two_norm() < tolerance;
    }

    /*!
     * \brief The imposed trace value at a position of the given cell of the trace shared
     *        with the given mortar domain: the cell's value for data per trace cell, the
     *        linear interpolation of the vertex values for data per trace vertex.
     */
    TraceValue traceValueAt(std::size_t mortarId, std::size_t traceCellIndex, const GlobalPosition& globalPos) const
    {
        const auto& data = traceData(mortarId);
        if (traceDataOrder_ == 0)
            return data[traceCellIndex];

        const auto& geometry = traceCellGeometry(mortarId, traceCellIndex);
        std::vector<Dune::FieldVector<Scalar, 1>> shapeValues;
        feCache_.get(geometry.type()).localBasis().evaluateFunction(geometry.local(globalPos), shapeValues);
        const auto& vertices = traceCellVertices_.at(mortarId)[traceCellIndex];
        TraceValue result(0.0);
        for (std::size_t i = 0; i < vertices.size(); ++i)
            result.axpy(shapeValues[i][0], data[vertices[i]]);
        return result;
    }

    //! Return the data imposed on the trace shared with the given mortar domain
    const TraceSolutionVector& traceData(std::size_t mortarId) const
    { return traceData_.at(mortarId); }

    //! Return the registered traces by the id of their mortar domain
    const std::unordered_map<std::size_t, std::shared_ptr<const Trace>>& traces() const
    { return traces_; }

    //! Return the entity couplings of the element with the given index
    const std::vector<EntityCouplingMap>& entityCouplingsOf(std::size_t eIdx) const
    { return eIdx < coupledEntities_.size() ? coupledEntities_[eIdx] : emptyEntityCouplings_; }

    //! Return the face couplings of the element with the given index
    const std::vector<ScvfCouplingMap>& faceCouplingsOf(std::size_t eIdx) const
    { return eIdx < coupledScvfs_.size() ? coupledScvfs_[eIdx] : emptyScvfCouplings_; }

    //! Return the boundary face couplings of the element with the given index
    const std::vector<BoundaryFaceCouplingMap>& boundaryFaceCouplingsOf(std::size_t eIdx) const
    { return eIdx < coupledBoundaryFaces_.size() ? coupledBoundaryFaces_[eIdx] : emptyBoundaryFaceCouplings_; }

    /*!
     * \brief Set whether the subdomain is solved without the external data of its own
     *        problem, which turns its solution into the action of the linear operator on
     *        the mortar data.
     * \note The problem keeps the responsibility for its own data; it asks with
     *       isHomogeneous(). Data this manager imposes is unaffected.
     */
    void setHomogeneous(bool value)
    { homogeneous_ = value; }

    //! Return true while the subdomain is solved without the external data of its problem
    bool isHomogeneous() const
    { return homogeneous_; }

    //! Return true if the given sub-control entity of the given element lies on a trace
    //! shared with a mortar domain
    template<typename SubControlEntity>
    bool isCoupled(const Element& element, const SubControlEntity& sce) const
    { return isOnMortarBoundary(*this, element, sce); }

    //! Return true if the given position lies on a trace cell shared with a mortar domain
    bool isCoupledAtPos(const GlobalPosition& globalPos) const
    {
        for (const auto& [id, trace] : traces_)
            if (!intersectingEntities(globalPos, trace->boundingBoxTree()).empty())
                return true;
        return false;
    }

    //! Return the mortar data imposed at the given sub-control entity of the given element,
    //! or at the given position on its boundary
    template<typename SubControlEntityOrPosition>
    decltype(auto) traceAt(const Element& element, const SubControlEntityOrPosition& sceOrPos) const
    { return traceVariablesAt(*this, element, sceOrPos); }

    //! Return the mortar data imposed at the given interpolation point of a boundary face
    template<typename FVElementGeometry, typename IpData>
    decltype(auto) traceAt(const FVElementGeometry& fvGeometry, const IpData& ipData) const
    {
        if constexpr (requires { ipData.scvfIndex(); })
            return traceVariablesAt(*this, fvGeometry.element(), fvGeometry.scvf(ipData.scvfIndex()));
        else
            return traceVariablesAt(*this, fvGeometry.element(), ipData.global());
    }

    /*!
     * \brief Return true if the degree of freedom of the given sub-control volume is the one
     *        to be pinned.
     *
     * A subdomain carrying mortar data on its entire boundary has no essential condition in
     * natural mode, so its problem is singular; fixing the first degree of freedom selects
     * the particular solution that vanishes there.
     */
    bool isPinnedDof(const Element& element, const SubControlVolume& scv) const
    { return couplingMode_ == CouplingMode::natural && scv.dofIndex() == 0 && isFloating(); }

private:
    using TraceCellVertices = Dune::ReservedVector<std::size_t, (1 << MortarGrid::dimension)>;
    using FECache = Dune::LagrangeLocalFiniteElementCache<typename GlobalPosition::value_type, Scalar, MortarGrid::dimension, 1>;

    std::shared_ptr<const GridGeometry> gridGeometry_;
    std::vector<std::vector<EntityCouplingMap>> coupledEntities_;
    std::vector<std::vector<ScvfCouplingMap>> coupledScvfs_;
    std::vector<std::vector<BoundaryFaceCouplingMap>> coupledBoundaryFaces_;
    std::unordered_map<std::size_t, std::vector<TraceCellGeometry>> traceCellGeometries_;
    std::unordered_map<std::size_t, std::vector<TraceCellVertices>> traceCellVertices_;
    std::unordered_map<std::size_t, std::shared_ptr<const Trace>> traces_;
    std::unordered_map<std::size_t, TraceSolutionVector> traceData_;
    FECache feCache_;
    std::size_t traceDataOrder_ = isCVFE ? 1 : 0;
    std::size_t dataVersion_ = 0;
    CouplingMode couplingMode_ = CouplingMode::essential;
    bool homogeneous_ = false;
    mutable std::optional<bool> isFloating_;

    static inline const std::vector<EntityCouplingMap> emptyEntityCouplings_{};
    static inline const std::vector<ScvfCouplingMap> emptyScvfCouplings_{};
    static inline const std::vector<BoundaryFaceCouplingMap> emptyBoundaryFaceCouplings_{};
};

#ifndef DOXYGEN
namespace Detail {

template<typename Traces, typename SubControlEntity>
std::size_t subControlEntityIndex(const SubControlEntity& sce)
{
    if constexpr (std::is_same_v<SubControlEntity, typename Traces::SubControlVolume>)
        return sce.dofIndex();
    else
        return sce.index();
}

//! The type of the sub-entities of the given kind of a subdomain grid geometry
template<typename GG, TraceEntity entity>
struct TraceSubEntity;

template<typename GG>
struct TraceSubEntity<GG, TraceEntity::subControlVolumeFace>
{ using type = typename GG::SubControlVolumeFace; };

template<typename GG>
struct TraceSubEntity<GG, TraceEntity::boundaryFace>
{ using type = typename GG::BoundaryFace; };

template<typename GG>
struct TraceSubEntity<GG, TraceEntity::subControlVolume>
{ using type = typename GG::SubControlVolume; };

} // end namespace Detail
#endif // DOXYGEN

/*!
 * \ingroup MortarCoupling
 * \brief Return true if the given sub-control entity of the given element lies on a trace
 *        shared with a mortar domain.
 */
template<typename GG, typename M, typename S, typename SubControlEntity>
bool isOnMortarBoundary(const FacetTraceCouplingManager<GG, M, S>& traces,
                        const typename FacetTraceCouplingManager<GG, M, S>::Element& element,
                        const SubControlEntity& sce)
{
    const auto eIdx = traces.gridGeometry().elementMapper().index(element);
    const auto sceIdx = Detail::subControlEntityIndex<FacetTraceCouplingManager<GG, M, S>>(sce);
    return std::ranges::any_of(traces.entityCouplingsOf(eIdx),
                               [&] (const auto& e) { return e.sceIndex == sceIdx; });
}

/*!
 * \ingroup MortarCoupling
 * \brief Return the imposed trace value at the integration point of the given boundary
 *        sub-control volume face: the value of the trace cell containing the face for data
 *        per trace cell, the trace function at the integration point for data per trace
 *        vertex.
 */
template<typename GG, typename M, typename S>
auto traceVariablesAt(const FacetTraceCouplingManager<GG, M, S>& traces,
                      const typename FacetTraceCouplingManager<GG, M, S>::Element& element,
                      const typename FacetTraceCouplingManager<GG, M, S>::SubControlVolumeFace& scvf)
{
    const auto eIdx = traces.gridGeometry().elementMapper().index(element);
    for (const auto& entry : traces.faceCouplingsOf(eIdx))
        if (entry.scvfIndex == scvf.index())
            return traces.traceValueAt(entry.mortarId, entry.traceCellIndex, scvf.ipGlobal());
    DUNE_THROW(Dune::InvalidStateException, "Given scvf does not overlap with a mortar domain.");
}

/*!
 * \ingroup MortarCoupling
 * \brief Return the imposed trace value associated with the given sub-control volume. For a
 *        control-volume finite element scheme this is the trace function at the degree of
 *        freedom: the mean over the trace vertices it coincides with for data per trace
 *        vertex, the mean over all trace cells containing it for data per trace cell, so
 *        that every element sharing the degree of freedom reads the same value. For a
 *        cell-centred scheme it is the mean over the coupled faces of the cell, which takes
 *        data per trace cell.
 */
template<typename GG, typename M, typename S>
auto traceVariablesAt(const FacetTraceCouplingManager<GG, M, S>& traces,
                      const typename FacetTraceCouplingManager<GG, M, S>::Element& element,
                      const typename FacetTraceCouplingManager<GG, M, S>::SubControlVolume& scv)
{
    using Traces = FacetTraceCouplingManager<GG, M, S>;
    const auto eIdx = traces.gridGeometry().elementMapper().index(element);
    typename Traces::TraceValue result(0);
    int count = 0;
    if constexpr (Traces::isCVFE)
    {
        if (traces.traceDataOrder() == 0)
        {
            for (const auto& [mortarId, trace] : traces.traces())
                for (const auto traceCellIndex : intersectingEntities(scv.dofPosition(), trace->boundingBoxTree()))
                {
                    result += traces.traceData(mortarId)[traceCellIndex];
                    ++count;
                }
        }
        else
            for (const auto& entry : traces.entityCouplingsOf(eIdx))
                if (entry.sceIndex == scv.dofIndex())
                {
                    result += traces.traceData(entry.mortarId)[entry.traceDofIndex];
                    ++count;
                }
    }
    else
    {
        if (traces.traceDataOrder() != 0)
            DUNE_THROW(Dune::NotImplemented, "Trace data per vertex read at the cell of a cell-centred scheme; read it at the faces");
        for (const auto& entry : traces.entityCouplingsOf(eIdx))
            if (entry.sceIndex == scv.dofIndex())
            {
                result += traces.traceData(entry.mortarId)[entry.traceDofIndex];
                ++count;
            }
    }
    if (count == 0)
        DUNE_THROW(Dune::InvalidStateException, "Given scv does not overlap with a mortar domain.");
    result /= count;
    return result;
}

/*!
 * \ingroup MortarCoupling
 * \brief Return the imposed trace value at the given position on the boundary of the given
 *        element: the trace function evaluated on the trace cell containing the position,
 *        independent of the relative resolution of trace and subdomain faces.
 */
template<typename GG, typename M, typename S>
auto traceVariablesAt(const FacetTraceCouplingManager<GG, M, S>& traces,
                      const typename FacetTraceCouplingManager<GG, M, S>::Element& element,
                      const typename FacetTraceCouplingManager<GG, M, S>::GlobalPosition& globalPos)
{
    using Traces = FacetTraceCouplingManager<GG, M, S>;
    using ctype = typename Traces::GlobalPosition::value_type;

    const auto eIdx = traces.gridGeometry().elementMapper().index(element);
    const typename Traces::BoundaryFaceCouplingMap* best = nullptr;
    bool bestContained = false;
    auto bestDistance = std::numeric_limits<ctype>::max();
    for (const auto& entry : traces.boundaryFaceCouplingsOf(eIdx))
    {
        const bool contained = traces.traceCellContains(entry.mortarId, entry.traceCellIndex, globalPos);
        const auto distance = (traces.traceCellGeometry(entry.mortarId, entry.traceCellIndex).center() - globalPos).two_norm();
        if (std::tie(bestContained, distance) < std::tie(contained, bestDistance))
        {
            best = &entry;
            bestContained = contained;
            bestDistance = distance;
        }
    }
    if (!best)
        DUNE_THROW(Dune::InvalidStateException, "Given element does not touch a mortar domain.");
    return traces.traceValueAt(best->mortarId, best->traceCellIndex, globalPos);
}

/*!
 * \ingroup MortarCoupling
 * \brief The mean of a quantity over each cell of the trace shared with the given mortar
 *        domain, assembled from the subdomain's boundary sub-control volume faces or its
 *        boundary faces.
 *
 * The field is bound per element and evaluated per coupled sub-entity of the chosen kind,
 * returning the integral of the quantity over that sub-entity; the sum over a trace cell's
 * sub-entities divided by the cell's measure is the coefficient of the quantity's trace in
 * the piecewise-constant space on the trace cells, which the mortar projectors consume.
 */
template<TraceEntity entity, typename GG, typename M, typename S, typename Field>
    requires Concept::FaceField<Field,
                                typename FacetTraceCouplingManager<GG, M, S>::Element,
                                typename Detail::TraceSubEntity<GG, entity>::type>
S traceValues(const FacetTraceCouplingManager<GG, M, S>& traces, std::size_t mortarId, const Field& field)
{
    static_assert(entity != TraceEntity::subControlVolume, "A facet trace is assembled from faces");
    using Extrusion = typename GG::Extrusion;
    S result(traces.numTraceCells(mortarId));
    result = 0.0;
    std::vector<typename S::field_type> measure(result.size(), 0.0);

    const auto couplingsOf = [&] (std::size_t eIdx) -> const auto& {
        if constexpr (entity == TraceEntity::subControlVolumeFace)
            return traces.faceCouplingsOf(eIdx);
        else
            return traces.boundaryFaceCouplingsOf(eIdx);
    };

    auto fvGeometry = localView(traces.gridGeometry());
    for (const auto& element : elements(traces.gridGeometry().gridView()))
    {
        const auto eIdx = traces.gridGeometry().elementMapper().index(element);
        if (std::ranges::none_of(couplingsOf(eIdx), [&] (const auto& e) { return e.mortarId == mortarId; }))
            continue;

        fvGeometry.bind(element);
        const auto bound = field.bind(element);
        const auto add = [&] (const auto& entry, const auto& subEntity) {
            result[entry.traceCellIndex] += bound(subEntity);
            measure[entry.traceCellIndex] += Extrusion::area(fvGeometry, subEntity);
        };
        for (const auto& entry : couplingsOf(eIdx))
        {
            if (entry.mortarId != mortarId)
                continue;
            if constexpr (entity == TraceEntity::subControlVolumeFace)
                add(entry, fvGeometry.scvf(entry.scvfIndex));
            else
                for (const auto& face : boundaryFaces(fvGeometry))
                    if (face.intersectionIndex() == entry.intersectionIndex)
                        add(entry, face);
        }
    }

    for (std::size_t i = 0; i < result.size(); ++i)
        if (measure[i] > 0.0)
            result[i] /= measure[i];
    return result;
}

/*!
 * \ingroup MortarCoupling
 * \brief The coupling manager of a subdomain with the given grid geometry on a mortar grid.
 *
 * A specialization defines the member type `type`, the coupling manager for the trace kind
 * the subdomain has on the mortar grid (see TraceTraits). The one given here serves traces
 * whose cells form a grid; a new trace kind specializes this trait next to its manager.
 */
template<typename GG, typename MortarGrid, typename MortarSolution>
struct SubDomainCouplingManagerTraits
{
    static_assert(HasTrace<GG, MortarGrid>,
                  "No mortar trace kind is defined for this subdomain grid geometry and mortar "
                  "grid. Include the header defining the trace kind of the subdomain.");
    static_assert(!HasTrace<GG, MortarGrid>,
                  "No subdomain coupling manager is defined for the trace kind of this subdomain. "
                  "Include the header defining the coupling manager of the trace kind.");
};

template<typename GG, typename MortarGrid, typename MortarSolution>
    requires Concept::GridMortarTrace<Trace<GG, MortarGrid>>
struct SubDomainCouplingManagerTraits<GG, MortarGrid, MortarSolution>
{ using type = FacetTraceCouplingManager<GG, MortarGrid, MortarSolution>; };

/*!
 * \ingroup MortarCoupling
 * \brief The coupling manager of a subdomain in a mortar-coupled model, selected by the trace
 *        kind of the subdomain. Subdomain problems compose this type.
 *
 * The selection is by the kind of the subdomain's trace, so it is decided once the trace
 * kind is known and does not depend on the discretization scheme of the subdomain. A
 * subdomain whose trace kind has no manager fails to compile here.
 */
template<typename GG, typename MortarGrid, typename MortarSolution>
using SubDomainCouplingManager = typename SubDomainCouplingManagerTraits<GG, MortarGrid, MortarSolution>::type;

} // end namespace Dumux::Mortar

#endif
