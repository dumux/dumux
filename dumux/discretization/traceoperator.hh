// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Classes to compute trace & normal traces of grid functions.
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
    { t.gridGeometryLocalView() };
};

//! Convenience implementation for trace operator contexts for finite-volume schemes
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

    const auto& gridGeometryLocalView() const { return fvGeometry_; }
    const auto& elemVolVars() const { return elemVolVars_; }
    const auto& elemFluxVarsCache() const { return elemFluxVarsCache_; }

 private:
    const SolutionVector& x_;
    typename GridVariables::GridGeometry::LocalView fvGeometry_;
    typename GridVariables::GridVolumeVariables::LocalView elemVolVars_;
    typename GridVariables::GridFluxVariablesCache::LocalView elemFluxVarsCache_;
};

#ifndef DOXYGEN
namespace TraceOperatorDetail {

template<typename GridGeometry>
struct DefaultContext
{
    typename GridGeometry::LocalView localView;
    void bind(const auto& element) { localView.bind(element); }
    const auto& gridGeometryLocalView() const { return localView; }
};

template<typename GridGeometry>
using VolumeType = std::remove_cvref_t<decltype(
    GridGeometry::Extrusion::area(
        std::declval<const typename GridGeometry::LocalView&>(),
        std::declval<const typename GridGeometry::SubControlVolumeFace&>()
    )
)>;

template<typename TraceGridGeometry, typename TraceGridMapper, typename TraceToDomainVertex>
std::vector<std::size_t> makeDomainToTraceVertexMap(const TraceGridGeometry& traceGridGeometry,
                                                    const TraceGridMapper& traceGridMapper,
                                                    const TraceToDomainVertex& traceToDomainVertex)
{
    std::vector<std::size_t> map;
    map.resize(traceGridMapper.domainGridGeometry().vertexMapper().size());
    for (const auto& v : vertices(traceGridGeometry.gridView()))
        map[
            traceGridMapper.domainGridGeometry().vertexMapper().index(traceToDomainVertex(v))
        ] = traceGridGeometry.vertexMapper().index(v);
    return map;
}

template<typename TraceGridGeometry, typename TraceGridMapper, typename LocalContext, typename F>
auto elementWiseAreaWeightedAverage(const TraceGridGeometry& traceGridGeometry,
                                    const TraceGridMapper& traceGridMapper,
                                    LocalContext context,
                                    F&& f)
{
    using GG = typename TraceGridMapper::DomainGridGeometry;
    using T = std::invoke_result_t<F, const typename GG::SubControlVolumeFace&, const LocalContext&>;

    Dune::BlockVector<T> result;
    result.resize(traceGridGeometry.elementMapper().size());
    for (const auto& traceElement : elements(traceGridGeometry.gridView()))
    {
        VolumeType<GG> area = 0.0;
        auto& entry = result[traceGridGeometry.elementMapper().index(traceElement)];
        for (const auto& domainElement : traceGridMapper.domainElementsAdjacentTo(traceElement) | std::views::take(1))
        {
            context.bind(domainElement);
            for (const auto& scvfIdx : traceGridMapper.domainScvfsAdjacentTo(traceElement, domainElement))
            {
                const auto& scvf = context.gridGeometryLocalView().scvf(scvfIdx);
                const auto scvfArea = GG::Extrusion::area(context.gridGeometryLocalView(), scvf);
                auto scvfContribution = f(scvf, context);
                scvfContribution *= scvfArea;
                entry += scvfContribution;
                area += scvfArea;
            }
        }
        entry /= area;
    }
    return result;
}

}  // namespace TraceOperatorDetail
#endif  // DOXYGEN

template<typename TraceGridGeometry,
         typename TraceGridMapper,
         typename Implementation>
class TraceOperatorBase
{
    using GridGeometry = typename TraceGridMapper::DomainGridGeometry;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<GridView::dimension>::Entity;

    using TraceGridView = typename TraceGridGeometry::GridView;
    using TraceVertex = TraceGridView::template Codim<TraceGridView::dimension>::Entity;

 public:
    using DefaultContext = TraceOperatorDetail::DefaultContext<GridGeometry>;

    TraceOperatorBase(std::shared_ptr<const TraceGridGeometry> traceGridGeometry,
                      std::shared_ptr<const TraceGridMapper> traceGridMapper)
    : traceGridGeometry_{std::move(traceGridGeometry)}
    , traceGridMapper_{std::move(traceGridMapper)}
    {}

    TraceOperatorBase(std::shared_ptr<const TraceGridGeometry> traceGridGeometry,
                      std::shared_ptr<const TraceGridMapper> traceGridMapper,
                      std::shared_ptr<const std::vector<std::size_t>> domainToTraceVertex)
    : traceGridGeometry_{std::move(traceGridGeometry)}
    , traceGridMapper_{std::move(traceGridMapper)}
    , domainToTraceVertex_{std::move(domainToTraceVertex)}
    {}

    template<std::invocable<const TraceVertex&> TraceToDomainVertex>
        requires(std::convertible_to<std::invoke_result_t<TraceToDomainVertex, const TraceVertex&>, const Vertex&>)
    TraceOperatorBase(std::shared_ptr<const TraceGridGeometry> traceGridGeometry,
                      std::shared_ptr<const TraceGridMapper> traceGridMapper,
                      const TraceToDomainVertex& mapper)
    : TraceOperatorBase{
        traceGridGeometry,
        traceGridMapper,
        std::make_shared<std::vector<std::size_t>>(
            TraceOperatorDetail::makeDomainToTraceVertexMap(*traceGridGeometry, *traceGridMapper, mapper)
        )
    }
    {}

    /*!
     * \brief Compute the trace coefficients for a given grid function.
     * \param f Functor to retrieve the grid function values at trace scvfs
     */
    template<std::invocable<const SubControlVolumeFace&, const DefaultContext&> F>
    auto apply(F&& f) const
    { return apply(std::forward<F>(f), DefaultContext{localView(traceGridMapper_->domainGridGeometry())}); }

    /*!
     * \brief Compute the trace coefficients for a given grid function.
     * \param f Functor to retrieve the grid function values at trace scvfs
     * \param context A context class that will be bound per visited domain element and passed to f
     * \note the context may be used to pre-bind elemVolVars or precompute other element-local variables.
     */
    template<typename F, TraceOperatorContext<Element> LocalContext>
        requires(std::is_invocable_v<F, const SubControlVolumeFace&, const LocalContext&>)
    auto apply(F&& f, LocalContext context) const
    { return static_cast<const Implementation&>(*this).apply_(std::forward<F>(f), std::move(context)); }

 protected:
    const auto& domainToTraceVertexMapPtr_() const
    { return domainToTraceVertex_; }

    const auto& domainToTraceVertexMap_() const
    {
        if (not domainToTraceVertex_)
            DUNE_THROW(Dune::InvalidStateException, "Trace operator has not been constructed with a trace vertex map");
        return *domainToTraceVertex_;
    }

    std::shared_ptr<const TraceGridGeometry> traceGridGeometry_;
    std::shared_ptr<const TraceGridMapper> traceGridMapper_;

 private:
    std::shared_ptr<const std::vector<std::size_t>> domainToTraceVertex_{nullptr};
};

//! Computes the trace of a grid function
template<typename TraceGridGeometry, typename TraceGridMapper>
class TraceOperator;

template<typename TGG, typename TGM, typename... Args>
TraceOperator(std::shared_ptr<TGG>, std::shared_ptr<TGM>, Args&&...)
-> TraceOperator<std::remove_cvref_t<TGG>, std::remove_cvref_t<TGM>>;

// Specialization for cell-centered schemes
template<typename TraceGridGeometry, typename TraceGridMapper>
    requires(
        TraceGridMapper::DomainGridGeometry::discMethod == DiscretizationMethods::cctpfa ||
        TraceGridMapper::DomainGridGeometry::discMethod == DiscretizationMethods::ccmpfa
    )
class TraceOperator<TraceGridGeometry, TraceGridMapper>
: public TraceOperatorBase<TraceGridGeometry, TraceGridMapper, TraceOperator<TraceGridGeometry, TraceGridMapper>>
{
    using Base = TraceOperatorBase<TraceGridGeometry, TraceGridMapper, TraceOperator<TraceGridGeometry, TraceGridMapper>>;

 public:
    using Base::Base;

 private:
    friend Base;

    template<typename F, typename LocalContext>
    auto apply_(F&& f, LocalContext context) const
    {
        return TraceOperatorDetail::elementWiseAreaWeightedAverage(
            *this->traceGridGeometry_,
            *this->traceGridMapper_,
            std::move(context),
            std::forward<F>(f)
        );
    }
};

// Specialization for the box scheme
template<typename TraceGridGeometry, typename TraceGridMapper>
    requires(TraceGridMapper::DomainGridGeometry::discMethod == DiscretizationMethods::box)
class TraceOperator<TraceGridGeometry, TraceGridMapper>
: public TraceOperatorBase<TraceGridGeometry, TraceGridMapper, TraceOperator<TraceGridGeometry, TraceGridMapper>>
{
    using Base = TraceOperatorBase<TraceGridGeometry, TraceGridMapper, TraceOperator<TraceGridGeometry, TraceGridMapper>>;
    using GG = TraceGridMapper::DomainGridGeometry;

 public:
    using Base::Base;

 private:
    friend Base;

    template<typename F, typename LocalContext>
    auto apply_(F&& f, LocalContext context) const
    {
        Dune::BlockVector<
            std::invoke_result_t<F, const typename GG::SubControlVolumeFace&, const LocalContext&>
        > result;
        result.resize(this->traceGridGeometry_->vertexMapper().size());
        for (const auto& traceElement : elements(this->traceGridGeometry_->gridView()))
        {
            for (const auto& domainElement : this->traceGridMapper_->domainElementsAdjacentTo(traceElement) | std::views::take(1))
            {
                context.bind(domainElement);
                for (const auto& scvfIdx : this->traceGridMapper_->domainScvfsAdjacentTo(traceElement, domainElement))
                {
                    const auto& scvf = context.gridGeometryLocalView().scvf(scvfIdx);
                    const auto& scv = context.gridGeometryLocalView().scv(scvf.insideScvIdx());
                    result[this->domainToTraceVertexMap_()[scv.dofIndex()]] = f(scvf, context);
                }
            }
        }
        return result;
    }
};


//! Computes the normal trace of a grid function
template<typename TraceGridGeometry, typename TraceGridMapper>
class NormalTraceOperator;

template<typename TGG, typename TGM, typename... Args>
NormalTraceOperator(std::shared_ptr<TGG>, std::shared_ptr<TGM>, Args&&...)
-> NormalTraceOperator<std::remove_cvref_t<TGG>, std::remove_cvref_t<TGM>>;

// Specialization for cell-centered schemes
template<typename TraceGridGeometry, typename TraceGridMapper>
    requires(
        TraceGridMapper::DomainGridGeometry::discMethod == DiscretizationMethods::cctpfa ||
        TraceGridMapper::DomainGridGeometry::discMethod == DiscretizationMethods::ccmpfa
    )
class NormalTraceOperator<TraceGridGeometry, TraceGridMapper>
: public TraceOperatorBase<TraceGridGeometry, TraceGridMapper, NormalTraceOperator<TraceGridGeometry, TraceGridMapper>>
{
    using Base = TraceOperatorBase<TraceGridGeometry, TraceGridMapper, NormalTraceOperator<TraceGridGeometry, TraceGridMapper>>;

 public:
    using Base::Base;

 private:
    friend Base;

    template<typename F, typename LocalContext>
    auto apply_(F&& f, LocalContext context) const
    {
        return TraceOperatorDetail::elementWiseAreaWeightedAverage(
            *this->traceGridGeometry_,
            *this->traceGridMapper_,
            std::move(context),
            std::forward<F>(f)
        );
    }
};

// Specialization for the box scheme
template<typename TraceGridGeometry, typename TraceGridMapper>
    requires(TraceGridMapper::DomainGridGeometry::discMethod == DiscretizationMethods::box)
class NormalTraceOperator<TraceGridGeometry, TraceGridMapper>
: public TraceOperatorBase<TraceGridGeometry, TraceGridMapper, NormalTraceOperator<TraceGridGeometry, TraceGridMapper>>
{
    using Base = TraceOperatorBase<TraceGridGeometry, TraceGridMapper, NormalTraceOperator<TraceGridGeometry, TraceGridMapper>>;

 public:
    using Base::Base;

 private:
    friend Base;

    template<typename F, typename LocalContext>
    auto apply_(F&& f, LocalContext context) const
    {
        return TraceOperatorDetail::elementWiseAreaWeightedAverage(
            *this->traceGridGeometry_,
            *this->traceGridMapper_,
            std::move(context),
            std::forward<F>(f)
        );
    }
};


//! Factory for different trace operators
template<typename TraceGridGeometry, typename TraceGridMapper>
class TraceOperatorFactory
: private TraceOperatorBase<TraceGridGeometry, TraceGridMapper, NormalTraceOperator<TraceGridGeometry, TraceGridMapper>>
{
    using Base = TraceOperatorBase<TraceGridGeometry, TraceGridMapper, NormalTraceOperator<TraceGridGeometry, TraceGridMapper>>;

 public:
    using Base::Base;

    //! Return the trace operator for this trace grid geometry & mapping
    TraceOperator<TraceGridGeometry, TraceGridMapper> traceOperator() const
    {
        auto ptr = this->domainToTraceVertexMapPtr_();
        if (ptr)
            return {this->traceGridGeometry_, this->traceGridMapper_, std::move(ptr)};
        return {this->traceGridGeometry_, this->traceGridMapper_};
    }

    //! Return the normal trace operator for this trace grid geometry & mapping
    NormalTraceOperator<TraceGridGeometry, TraceGridMapper> normalTraceOperator() const
    {
        auto ptr = this->domainToTraceVertexMapPtr_();
        if (ptr)
            return {this->traceGridGeometry_, this->traceGridMapper_, std::move(ptr)};
        return {this->traceGridGeometry_, this->traceGridMapper_};
    }
};

template<typename TGG, typename TGM, typename... Args>
TraceOperatorFactory(std::shared_ptr<TGG>, std::shared_ptr<TGM>, Args&&...)
-> TraceOperatorFactory<std::remove_cvref_t<TGG>, std::remove_cvref_t<TGM>>;

} // end namespace Dumux

#endif
