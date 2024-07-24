// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MultiDomain ???
 * \brief Extract values on the domain boundary from finite-volume discretizations.
 * \todo This header should probably go to dumux/discretization?
 */
#ifndef DUMUX_MULTIDOMAIN_TRACE_HH
#define DUMUX_MULTIDOMAIN_TRACE_HH

#include <vector>
#include <memory>
#include <utility>
#include <concepts>
#include <type_traits>
#include <algorithm>
#include <unordered_map>
#include <ranges>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dumux/common/typetraits/problem.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/geometry/geometryintersection.hh>

namespace Dumux {

#ifndef DOXYGEN
namespace TraceDetail {

    template<typename Geo1, typename Geo2>
    bool intersect(const Geo1& geo1, const Geo2& geo2)
    {
        using Algo = GeometryIntersection<Geo1, Geo2>;
        typename Algo::Intersection r;
        return Algo::intersection(geo1, geo2, r);
    }

    template<typename GridVariables>
    using DefaultTraceGrid = Dune::FoamGrid<
        GridVariables::GridGeometry::GridView::dimension - 1,
        GridVariables::GridGeometry::GridView::dimensionworld
    >;

}  // namespace TraceDetail
#endif  // DOXYGEN


/*!
 * \ingroup MultiDomain ???
 * \brief Extract values on the domain boundary from finite-volume discretizations.
 */
template<typename GridVariables, typename TraceGrid = TraceDetail::DefaultTraceGrid<GridVariables>>
class FVTrace
{
    using GridGeometry = typename GridVariables::GridGeometry;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;

    using GridView = typename GridGeometry::GridView;
    using GridIntersection = typename GridView::Intersection;
    using GridElement = typename GridView::template Codim<0>::Entity;
    using Coordinate = typename GridElement::Geometry::GlobalCoordinate;

    using TraceElementToScvfElementIndices = std::unordered_map<std::size_t, std::vector<std::size_t>>;

 public:
    using TraceGridView = typename TraceGrid::LeafGridView;
    using TraceEntityMapper = Dune::MultipleCodimMultipleGeomTypeMapper<TraceGridView>;

    template<typename TraceIndicator>
    explicit FVTrace(std::shared_ptr<const GridVariables> gv, TraceIndicator&& indicator)
    : domainGridVariables_{gv}
    { update(std::forward<TraceIndicator>(indicator)); }

    //! Recompute the trace grid and mappings (e.g. after grid adaptation or with different indicator)
    template<std::invocable<const GridIntersection&> TraceIndicator>
        requires(std::is_convertible_v<std::invoke_result_t<TraceIndicator, const GridIntersection&>, bool>)
    void update(TraceIndicator&& indicator)
    {
        traceGridFactory_ = std::make_unique<Dune::GridFactory<TraceGrid>>();

        std::unordered_map<std::size_t, std::size_t> insertedVertices;
        std::vector<unsigned int> localCornerStorage;
        std::size_t insertedElementCount = 0;

        traceToScvfIndices_.resize(domainGridView_().size(0));
        auto fvGeometry = localView(domainGridGeometry_());
        std::ranges::for_each(elements(domainGridView_()), [&] (const auto& element) {
            const auto eIdx = domainGridGeometry_().elementMapper().index(element);
            const auto& refElement = Dune::referenceElement(element);
            const auto& elemGeo = element.geometry();

            auto getFVGeometry = [&, isBound=false] () mutable -> const typename GridGeometry::LocalView& {
                if (!isBound) { fvGeometry.bindElement(element); isBound = true; }
                return fvGeometry;
            };

            std::ranges::for_each(
                intersections(domainGridView_(), element)
                | std::views::filter([&] (const auto& is) { return is.boundary(); })
                | std::views::filter([&] (const auto& is) { return indicator(is); }),
                [&] (const auto& intersection) {
                    const auto fvGeometry = getFVGeometry();
                    const auto& isGeo = intersection.geometry();

                    localCornerStorage.clear();
                    localCornerStorage.reserve(isGeo.corners());
                    std::ranges::for_each(std::views::iota(0, isGeo.corners()), [&] (int c) {
                        const auto vIdxLocal = refElement.subEntity(intersection.indexInInside(), 1, c, GridView::dimension);
                        const auto vIdxGlobal = domainGridGeometry_().vertexMapper().subIndex(element, vIdxLocal, GridView::dimension);
                        localCornerStorage.push_back([&] () {
                            auto it = insertedVertices.find(vIdxGlobal);
                            if (it == insertedVertices.end())
                            {
                                const auto currentVertexIndex = insertedVertices.size();
                                traceGridFactory_->insertVertex(elemGeo.global(refElement.position(vIdxLocal, GridView::dimension)));
                                it = insertedVertices.emplace(static_cast<std::size_t>(vIdxGlobal), currentVertexIndex).first;
                            }
                            return it->second;
                        } ());
                    });

                    traceGridFactory_->insertElement(isGeo.type(), localCornerStorage);
                    insertedElementCount++;
                    std::ranges::for_each(scvfs(fvGeometry), [&] (const auto& scvf) {
                        if (TraceDetail::intersect(fvGeometry.geometry(scvf), isGeo))
                            traceToScvfIndices_[eIdx][insertedElementCount - 1].push_back(scvf.index());
                    });
            });
        });

        traceGrid_ = traceGridFactory_->createGrid();
        traceElementMapper_ = std::make_unique<TraceEntityMapper>(gridView(), Dune::mcmgElementLayout());
        traceVertexMapper_ = std::make_unique<TraceEntityMapper>(gridView(), Dune::mcmgVertexLayout());

        // map insertion to mapper indices for scvf map
        std::vector<std::size_t> insertionToActualIndex(traceGrid_->leafGridView().size(0));
        for (const auto& e : elements(traceGrid_->leafGridView()))
            insertionToActualIndex[traceGridFactory_->insertionIndex(e)] = traceElementMapper_->index(e);
        for (auto& map : traceToScvfIndices_) {
            TraceElementToScvfElementIndices new_map;
            for (const auto& [traceElemInsertionIndex, scvfIndices] : map)
                new_map[insertionToActualIndex[traceElemInsertionIndex]] = scvfIndices;
            map = std::move(new_map);
        }
    }

    //! Return the grid view representing the trace
    TraceGridView gridView() const
    { return traceGrid_->leafGridView(); }

    //! Return the index mapper for trace grid elements
    const TraceEntityMapper& elementMapper() const
    { return *traceElementMapper_; }

    //! Return the index mapper for trace grid vertices
    const TraceEntityMapper& vertexMapper() const
    { return *traceVertexMapper_; }

    //! Return the index of the trace element that the given scvf coincides with
    std::size_t traceElement(const GridElement& element, const SubControlVolumeFace& scvf) const
    {
        const auto eIdx = domainGridGeometry_().elementMapper().index(element);
        for (const auto& [traceElemIdx, scvfIndices] : traceToScvfIndices_.at(eIdx))
            if (scvfIndices.contains(scvf.index()))
                return traceElemIdx;
        DUNE_THROW(Dune::InvalidStateException, "Could not find trace element for the given scvf.");
    }

    //! Assemble variables on the trace
    template<typename SolutionVector, typename FaceVariableFunctor>
        requires(std::invocable<FaceVariableFunctor,
            const GridElement&,
            const typename GridGeometry::LocalView&,
            const typename GridVariables::GridVolumeVariables::LocalView&,
            const typename GridVariables::GridFluxVariablesCache::LocalView&,
            const SubControlVolumeFace&
        >)
    SolutionVector assemble(const SolutionVector& x, FaceVariableFunctor&& vars) const {
        SolutionVector trace;
        trace.resize(gridView().size(0));
        trace = 0;
        for (const auto& e : elements(domainGridView_())) {
            const auto eIdx = domainGridGeometry_().elementMapper().index(e);
            const auto& map = traceToScvfIndices_.at(eIdx);
            if (map.empty())
                continue;

            const auto fvGeometry = localView(domainGridGeometry_()).bind(e);
            const auto elemVolVars = localView(domainGridVariables_->curGridVolVars()).bind(e, fvGeometry, x);
            const auto elemFluxVarsCache = localView(domainGridVariables_->gridFluxVarsCache()).bind(e, fvGeometry, elemVolVars);
            for (const auto& [traceElemIdx, scvfIndices] : map) {
                if (!scvfIndices.empty()) {
                    typename GridVariables::Scalar totalArea = 0;
                    for (const auto scvfIdx : scvfIndices) {
                        const auto& scvf = fvGeometry.scvf(scvfIdx);
                        const auto area = Extrusion_t<GridGeometry>::area(fvGeometry, scvf);
                        trace[traceElemIdx] += area*vars(e, fvGeometry, elemVolVars, elemFluxVarsCache, scvf);
                        totalArea += area;
                    }
                    if (!scvfIndices.empty())
                        trace[traceElemIdx] /= totalArea;
                }
            }
        }
        return trace;
    }

 private:
    const auto& domainGridGeometry_() const { return domainGridVariables_->gridGeometry(); }
    auto domainGridView_() const { return domainGridGeometry_().gridView(); }

    std::shared_ptr<const GridVariables> domainGridVariables_;

    std::unique_ptr<Dune::GridFactory<TraceGrid>> traceGridFactory_;
    std::unique_ptr<TraceGrid> traceGrid_;

    std::unique_ptr<TraceEntityMapper> traceElementMapper_;
    std::unique_ptr<TraceEntityMapper> traceVertexMapper_;

    std::vector<TraceElementToScvfElementIndices> traceToScvfIndices_;
};

//! CTAD guide
template<typename GridVariables, typename Indicator>
FVTrace(std::shared_ptr<GridVariables>, Indicator&&) -> FVTrace<std::remove_const_t<GridVariables>>;


//! Trace operator (combines a trace and a given assembler function)
template<typename Trace, typename FaceAssembler>
class TraceOperator
{
 public:
    using TraceGridView = typename Trace::TraceGridView;
    using TraceEntityMapper = typename Trace::TraceEntityMapper;

    TraceOperator(Trace&& trace, FaceAssembler&& assembler)
    : trace_{std::move(trace)}
    , assembler_{std::move(assembler)}
    {}

    //! Return the grid view representing the trace
    TraceGridView gridView() const
    { return trace_.leafGridView(); }

    //! Return the index mapper for trace grid elements
    const TraceEntityMapper& elementMapper() const
    { return trace_.elementMapper(); }

    //! Return the index mapper for trace grid vertices
    const TraceEntityMapper& vertexMapper() const
    { return trace_.vertexMapper(); }

    //! Assemble variables on the trace according to the stored assembler
    template<typename SolutionVector>
    SolutionVector assemble(const SolutionVector& x) const {
        return trace_.assemble(x, assembler_);
    }

 private:
    Trace trace_;
    FaceAssembler assembler_;
};

//! Return the default advective flux function for trace flux assembly in cellcentered schemes (TODO: name should be trace related)
template<typename FluxVariables, typename Problem>
inline constexpr auto defaultCCAdvectiveFluxFunction(const Problem& problem)
{
    using GG = typename ProblemTraits<Problem>::GridGeometry;
    return [&] (
        const auto& element,
        const auto& fvGeometry,
        const auto& elemVolVars,
        const auto& elemFluxVarsCache,
        const auto& scvf
    ) {
        using NumEqVector = std::remove_cvref_t<decltype(problem.neumann(element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf))>;
        if (scvf.boundary() && problem.boundaryTypes(element, scvf).hasOnlyNeumann())
        {
            auto flux = problem.neumann(element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf);
            flux *= elemVolVars[fvGeometry.scv(scvf.insideScvIdx())].extrusionFactor();
            return flux;
        }

        NumEqVector result;
        FluxVariables fluxVars;
        fluxVars.init(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
        for (int i = 0; i < result.size(); ++i)
            result[i] = fluxVars.advectiveFlux(i, [&] (const auto& vv) { return vv.mobility(i); })
                        /Extrusion_t<GG>::area(fvGeometry, scvf);
        return result;
    };
}

//! TODO: This should go somewhere 1p-related
template<typename FluxVariables, typename Problem>
inline constexpr auto tracePressureFunctionTpfaOneP(const Problem& problem)
{
    using GG = typename ProblemTraits<Problem>::GridGeometry;
    return [&] (
        const auto& element,
        const auto& fvGeometry,
        const auto& elemVolVars,
        const auto& elemFluxVarsCache,
        const auto& scvf
    ) {
        const auto insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto fromFlux = [&] (auto flux) {
            const auto ti = elemFluxVarsCache[scvf].advectionTij()
                            *insideVolVars.density()
                            *insideVolVars.mobility();
            flux -= ti*insideVolVars.pressure();
            return -1.0*flux/ti;
        };

        if (scvf.boundary())
        {
            if (problem.boundaryTypes(element, scvf).hasOnlyNeumann())
                return fromFlux(
                    problem.neumann(element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf)[0]
                    *Extrusion_t<GG>::area(fvGeometry, scvf)
                    *insideVolVars.extrusionFactor()
                );
            else
                return elemVolVars[scvf.outsideScvIdx()].pressure();
        }

        FluxVariables fluxVars;
        fluxVars.init(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
        return fromFlux(fluxVars.advectiveFlux(0, [&] (const auto& vv) { return vv.mobility(0); }));
    };
}

} // end namespace Dumux

#endif
