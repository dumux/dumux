// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Class to assemble primary variables and/or fluxes on boundaries.
 */
#ifndef DUMUX_DISCRETIZATION_TRACE_OPERATOR_HH
#define DUMUX_DISCRETIZATION_TRACE_OPERATOR_HH

#include <config.h>

#include <memory>
#include <type_traits>

#include <dune/istl/bvector.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/facetgridmapper.hh>

namespace Dumux {

//! Concept for bindable local contexts (TODO: move this?)
template<typename T, typename A>
concept TraceOperatorContext = requires(T& t, const A& arg) {
    { t.bind(arg) };
    { t.fvGeometry() };
};

//! Default implementation for trace operator contexts for finite-volume schemes
template<typename GridVariables, typename SolutionVector>
class FVTraceOperatorContext
{
 public:
    FVTraceOperatorContext(const GridVariables& gv, const SolutionVector& x)
    : x_{x}
    , fvGeometry_{localView(gv.curGridVolVars().problem().gridGeometry())}
    , elemVolVars_{localView(gv.curGridVolVars())}
    , elemFluxVarsCache_{localView(gv.gridFluxVarsCache())}
    {}

    void bind(const auto& element)
    {
        fvGeometry_.bind(element);
        elemVolVars_.bind(element, fvGeometry_, x_);
        elemFluxVarsCache_.bind(element, fvGeometry_, elemVolVars_);
    }

    const auto& fvGeometry() const { return fvGeometry_; }
    const auto& elemVolVars() const { return elemVolVars_; }
    const auto& elemFluxVarsCache() const { return elemFluxVarsCache_; }

 private:
    const SolutionVector& x_;
    typename GridVariables::GridGeometry::LocalView fvGeometry_;
    typename GridVariables::GridVolumeVariables::LocalView elemVolVars_;
    typename GridVariables::GridFluxVariablesCache::LocalView elemFluxVarsCache_;
};

/*!
 * \ingroup Discretization
 * \brief Assembles primary variables and/or fluxes on traces.
 * \note This assumes that the provided trace grid only contains boundary facets. If it contains interior
 *       ones, continuity of variables/fluxes is assumed, and only one side of the facets will be visited.
 */
template<typename TraceGridGeometry, typename GridGeometry>
class FVTraceOperator {
    using TraceGridView = typename TraceGridGeometry::GridView;
    static constexpr int dim = GridGeometry::GridView::dimension;
    static constexpr int traceDim = TraceGridGeometry::GridView::dimension;

    struct DefaultContext
    {
        typename GridGeometry::LocalView localView;

        void bind(const auto& element) { localView.bind(element); }
        const auto& fvGeometry() const { return localView; }
    };

 public:
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using Vertex = typename GridGeometry::GridView::template Codim<dim>::Entity;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using TraceVertex = typename TraceGridGeometry::GridView::template Codim<traceDim>::Entity;

    template<std::invocable<const TraceVertex&> TraceToDomainVertexMap>
        requires(std::convertible_to<std::invoke_result_t<TraceToDomainVertexMap, const TraceVertex&>, const Vertex&>)
    FVTraceOperator(std::shared_ptr<const TraceGridGeometry> traceGridGeometry,
                    std::shared_ptr<const FVFacetGridMapper<TraceGridView, GridGeometry>> traceGridMapper,
                    const TraceToDomainVertexMap& traceToDomainVertex)
    : FVTraceOperator{std::move(traceGridGeometry), std::move(traceGridMapper)}
    {
        domainToTraceVertex_.resize(traceGridMapper_->domainGridGeometry().vertexMapper().size());
        for (const auto& v : vertices(traceGridGeometry_->gridView()))
            domainToTraceVertex_[
                traceGridMapper_->domainGridGeometry().vertexMapper().index(traceToDomainVertex(v))
            ] = traceGridGeometry_->vertexMapper().index(v);
    }

    FVTraceOperator(std::shared_ptr<const TraceGridGeometry> traceGridGeometry,
                    std::shared_ptr<const FVFacetGridMapper<TraceGridView, GridGeometry>> traceGridMapper)
    : traceGridGeometry_{std::move(traceGridGeometry)}
    , traceGridMapper_{std::move(traceGridMapper)}
    {}

    /*!
     * \brief Assemble vertex-wise variables on the trace from overlapping sub-control volumes (CVFE methods only).
     * \param f Functor to assemble the desired variables in boundary scvs
     */
    template<std::invocable<const SubControlVolume&> F>
    auto assembleScvVariables(F&& f) const
    {
        return assembleScvVariables(
            [&] (const auto& scv, const auto&) { return f(scv); },
            DefaultContext{localView(traceGridMapper_->domainGridGeometry())}
        );
    }

    /*!
     * \brief Assemble vertex-wise variables on the trace from overlapping sub-control volumes (CVFE methods only).
     * \param f Functor to assemble the desired variables at boundary scvs
     * \param context A context class that will be bound per visited domain element and passed to f
     * \note the context may be used to pre-bind elemVolVars or precompute other element-local variables.
     */
    template<typename F, TraceOperatorContext<Element> LocalContext>
        requires(
            std::is_invocable_v<F, const SubControlVolume&, const LocalContext&> and
            DiscretizationMethods::isCVFE<typename GridGeometry::DiscretizationMethod>
        )
    auto assembleScvVariables(F&& f, LocalContext context) const
    {
        static_assert(
            DiscretizationMethods::isCVFE<typename GridGeometry::DiscretizationMethod>,
            "Scv-wise assembly is only available with CVFE methods"
        );

        if (domainToTraceVertex_.empty())
            DUNE_THROW(
                Dune::InvalidStateException,
                "Trace operator has not been constructed with a provided vertex mapping. This is required for scv-wise trace assembly."
            );

        using T = std::invoke_result_t<F, const SubControlVolume&, const LocalContext&>;
        Dune::BlockVector<T> result;
        result.resize(traceGridGeometry_->vertexMapper().size());
        for (const auto& traceElement : elements(traceGridGeometry_->gridView()))
            for (const auto& domainElement : traceGridMapper_->domainElementsAdjacentTo(traceElement) | std::views::take(1))
            {
                context.bind(domainElement);
                for (const auto& scvfIdx : traceGridMapper_->domainScvfsAdjacentTo(traceElement, domainElement))
                {
                    const auto& scv = context.fvGeometry().scv(context.fvGeometry().scvf(scvfIdx).insideScvIdx());
                    result[domainToTraceVertex_[scv.dofIndex()]] = f(scv, context);
                }
            }
        return result;
    }

    /*!
     * \brief Assemble element-wise variables on the trace from overlapping sub-control volume faces
     * \param f Functor to assemble the desired variables at boundary scvfs
     * \param context A context class that will be bound per visited domain element and passed to f
     */
    template<std::invocable<const SubControlVolumeFace&> F>
    auto assembleScvfVariables(F&& f) const
    {
        return assembleScvfVariables(
            [&] (const auto& scvf, const auto&) { return f(scvf); },
            DefaultContext{localView(traceGridMapper_->domainGridGeometry())}
        );
    }

    /*!
     * \brief Assemble element-wise variables on the trace from overlapping sub-control volume faces
     * \param f Functor to assemble the desired variables at boundary scvfs
     * \param context A context class that will be bound per visited domain element and passed to f
     * \note the context may be used to pre-bind elemVolVars or precompute other element-local variables.
     */
    template<typename F, TraceOperatorContext<Element> LocalContext>
        requires(std::is_invocable_v<F, const SubControlVolumeFace&, const LocalContext&>)
    auto assembleScvfVariables(F&& f, LocalContext context) const
    {
        using T = std::invoke_result_t<F, const SubControlVolumeFace&, const LocalContext&>;

        Dune::BlockVector<T> result;
        result.resize(traceGridGeometry_->elementMapper().size());
        for (const auto& traceElement : elements(traceGridGeometry_->gridView()))
        {
            const auto traceElementIndex = traceGridGeometry_->elementMapper().index(traceElement);

            result[traceElementIndex] = 0;
            double area = 0.0;
            for (const auto& domainElement : traceGridMapper_->domainElementsAdjacentTo(traceElement) | std::views::take(1))
            {
                context.bind(domainElement);
                for (const auto& scvfIdx : traceGridMapper_->domainScvfsAdjacentTo(traceElement, domainElement))
                {
                    const auto& scvf = context.fvGeometry().scvf(scvfIdx);
                    const auto scvfArea = GridGeometry::Extrusion::area(context.fvGeometry(), scvf);
                    T scvfContribution = f(scvf, context);
                    scvfContribution *= scvfArea;
                    result[traceElementIndex] += scvfContribution;
                    area += scvfArea;
                }
            }
            result[traceElementIndex] /= area;
        }
        return result;
    }

 private:
    std::shared_ptr<const TraceGridGeometry> traceGridGeometry_;
    std::shared_ptr<const FVFacetGridMapper<TraceGridView, GridGeometry>> traceGridMapper_;
    std::vector<std::size_t> domainToTraceVertex_;
};

template<typename TGG, typename Mapper, typename... Args>
FVTraceOperator(std::shared_ptr<TGG>, std::shared_ptr<Mapper>, Args&&...)
-> FVTraceOperator<std::remove_const_t<TGG>, typename std::remove_const_t<Mapper>::DomainGridGeometry>;

} // end namespace Dumux

#endif
