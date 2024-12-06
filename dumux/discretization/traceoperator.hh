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
#include <dumux/discretization/facetgrid.hh>

namespace Dumux {

//! Concept for bindable local contexts
template<typename T, typename A>
concept BindableWith = requires(T& t, const A& arg) {
    { t.bind(arg) };
};

// TODO: Default context for grid variables
template<typename GridVariables, typename SolutionVector>
class FVContext
{
 public:
    FVContext(const GridVariables& gv, const SolutionVector& x)
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
 * \note This assumes that the provided FVFacetGrid only contains boundary facets. If it contains interior
 *       ones, continuity of variables/fluxes is assumed, and only one side of the facets will be visited.
 */
template<typename BoundaryGrid, typename GridGeometry>
class FVTraceOperator {
    struct DefaultContext
    {
        typename GridGeometry::LocalView localView;

        void bind(const auto& element) { localView.bind(element); }
        const auto& fvGeometry() const { return localView; }
    };

 public:
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;

    FVTraceOperator(std::shared_ptr<const FVFacetGrid<BoundaryGrid, GridGeometry>> traceGrid)
    : traceGrid_{std::move(traceGrid)}
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
            DefaultContext{localView(*domainGridGeometry())}
        );
    }

    /*!
     * \brief Assemble vertex-wise variables on the trace from overlapping sub-control volumes (CVFE methods only).
     * \param f Functor to assemble the desired variables at boundary scvs
     * \param context A context class that will be bound per visited domain element and passed to f
     * \note the context may be used to pre-bind elemVolVars or precompute other element-local variables.
     */
    template<typename F, BindableWith<Element> LocalContext>
        requires(
            std::is_invocable_v<F, const SubControlVolume&, const LocalContext&> and
            DiscretizationMethods::isCVFE<typename GridGeometry::DiscretizationMethod>
        )
    auto assembleScvVariables(F&& f, LocalContext context) const
    {
        static_assert(
            requires(const LocalContext& c) {
                { c.fvGeometry() } -> std::convertible_to<const typename GridGeometry::LocalView&>;
            },
            "LocalContext is required to expose a local grid geometry view via context.fvGeometry()"
        );

        std::vector<std::size_t> domainToFacetVertexIndex(domainGridGeometry()->vertexMapper().size(), 0);
        for (const auto& v : vertices(traceGridView()))
            domainToFacetVertexIndex[traceGrid_->domainVertexIndexOf(v)] = traceGrid_->vertexMapper().index(v);

        using T = std::invoke_result_t<F, const SubControlVolume&, const LocalContext&>;
        Dune::BlockVector<T> result;
        result.resize(traceGrid_->vertexMapper().size());
        for (const auto& traceElement : elements(traceGrid_->gridView()))
            for (const auto& domainElement : traceGrid_->domainElementsAdjacentTo(traceElement) | std::views::take(1))
            {
                context.bind(domainElement);
                for (const auto& scvfIdx : traceGrid_->domainScvfsAdjacentTo(traceElement, domainElement))
                {
                    const auto& scv = context.fvGeometry().scv(context.fvGeometry().scvf(scvfIdx).insideScvIdx());
                    result[domainToFacetVertexIndex[scv.dofIndex()]] = f(scv, context);
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
            DefaultContext{localView(*domainGridGeometry())}
        );
    }

    /*!
     * \brief Assemble element-wise variables on the trace from overlapping sub-control volume faces
     * \param f Functor to assemble the desired variables at boundary scvfs
     * \param context A context class that will be bound per visited domain element and passed to f
     * \note the context may be used to pre-bind elemVolVars or precompute other element-local variables.
     */
    template<typename F, BindableWith<Element> LocalContext>
        requires(std::is_invocable_v<F, const SubControlVolumeFace&, const LocalContext&>)
    auto assembleScvfVariables(F&& f, LocalContext context) const
    {
        static_assert(
            requires(const LocalContext& c) {
                { c.fvGeometry() } -> std::convertible_to<const typename GridGeometry::LocalView&>;
            },
            "LocalContext is required to expose a local grid geometry view via context.fvGeometry()"
        );

        using T = std::invoke_result_t<F, const SubControlVolumeFace&, const LocalContext&>;

        Dune::BlockVector<T> result;
        result.resize(traceGrid_->elementMapper().size());
        for (const auto& traceElement : elements(traceGridView()))
        {
            const auto traceElementIndex = traceGrid_->elementMapper().index(traceElement);

            result[traceElementIndex] = 0;
            double area = 0.0;
            for (const auto& domainElement : traceGrid_->domainElementsAdjacentTo(traceElement) | std::views::take(1))
            {
                context.bind(domainElement);
                for (const auto& scvfIdx : traceGrid_->domainScvfsAdjacentTo(traceElement, domainElement))
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

    auto traceGridView() const
    { return traceGrid_->gridView(); }

    std::shared_ptr<const GridGeometry> domainGridGeometry() const
    { return traceGrid_->domainGridGeometry(); }

 private:
    std::shared_ptr<const FVFacetGrid<BoundaryGrid, GridGeometry>> traceGrid_;
};

template<typename BG, typename GG>
FVTraceOperator(std::shared_ptr<FVFacetGrid<BG, GG>>) -> FVTraceOperator<BG, GG>;

} // end namespace Dumux

#endif
