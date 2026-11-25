// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Core
 * \brief Point source types and helpers for handling point sources.
 */
#ifndef DUMUX_POINTSOURCES_HH
#define DUMUX_POINTSOURCES_HH

#include <algorithm>
#include <map>
#include <ranges>
#include <span>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/reservedvector.hh>

#include <dumux/common/numeqvector.hh>
#include <dumux/discretization/localview.hh>
#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/geometry/intersectingentities.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/concepts.hh>
#include <dumux/discretization/elementscontainingpoint.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup Core
 * \brief Empty point sources.
 *
 * Provides the full `contexts`/`eval` interface expected by assembly,
 * but `contexts` always returns an empty span so nothing is ever evaluated.
 */
struct EmptyPointSources
{
    //! Dummy context type satisfying the assembly interface.
    struct Context{};

    constexpr bool empty() const { return true; }

    template<class ElementDiscretization>
    std::span<const Context> contexts(const ElementDiscretization&,
                                      const typename ElementDiscretization::SubControlVolume&) const
    { return {}; }

    template<class ElementDiscretization>
    std::span<const Context> contexts(const ElementDiscretization&,
                                      const typename ElementDiscretization::Element&) const
    { return {}; }

    template<class ElementDiscretization, class ElementVariables>
    auto eval(const ElementDiscretization&, const ElementVariables&, const Context&) const
    {
        using NumEqVector = Dumux::NumEqVector<typename ElementVariables::Variables::PrimaryVariables>;
        return NumEqVector{};
    }
};

/*!
 * \ingroup Core
 * \brief Identifies a point source by its position in space.
 *        The embeddings value is set during map construction to
 *        account for sources shared between several geometrical entities (scvs).
 * \tparam PositionType the position type
 */
template<class PositionType>
class PointSourceContext
{
public:
    using GlobalPosition = PositionType;

    PointSourceContext() = default;

    PointSourceContext(GlobalPosition position, std::size_t embeddings = 1)
    : position_(std::move(position))
    , embeddings_(embeddings)
    {}

    const GlobalPosition& position() const
    { return position_; }

    std::size_t embeddings() const
    { return embeddings_; }

    void setEmbeddings(std::size_t embeddings)
    { embeddings_ = embeddings; }

private:
    GlobalPosition position_;
    std::size_t embeddings_ = 1;
};

/*!
 * \ingroup Core
 * \brief Point source context enriched with a unique identifier.
 * \tparam Context the base context type
 * \tparam IdType the source identifier type
 */
template<class Context, class IdType>
class IdPointSourceContext : public Context
{
public:
    using Id = IdType;

    IdPointSourceContext() = default;
    IdPointSourceContext(const Context& ctx, IdType id)
    : Context(ctx)
    , id_(id)
    {}

    IdType id() const
    { return id_; }

private:
    IdType id_{};
};

/*!
 * \ingroup Core
 * \brief Point source context enriched with the sub-control volume index.
 */
template<class Context>
class ScvPointSourceContext : public Context
{
public:
    ScvPointSourceContext() = default;
    ScvPointSourceContext(const Context& ctx, std::size_t scvIdx)
    : Context(ctx)
    , scvIndex_(scvIdx)
    {}

    std::size_t scvIndex() const
    { return scvIndex_; }

private:
    std::size_t scvIndex_ = 0;
};

/*!
 * \ingroup Core
 * \brief Point source context enriched with the index of the element that contains the point.
 *
 * Contexts of this kind are located by `makePointSourceMap` without a search of the grid.
 * A point shared by several elements is given as one context per element,
 * each with the number of those elements as embeddings.
 */
template<class Context>
class ElementIndexPointSourceContext : public Context
{
public:
    ElementIndexPointSourceContext() = default;
    ElementIndexPointSourceContext(const Context& ctx, std::size_t eIdx)
    : Context(ctx)
    , elementIndex_(eIdx)
    {}

    std::size_t elementIndex() const
    { return elementIndex_; }

private:
    std::size_t elementIndex_ = 0;
};

/*!
 * \ingroup Core
 * \brief Point source context enriched with the position in the local coordinates of an element.
 * \tparam Context the base context type
 * \tparam LocalPositionType the position type in local element coordinates
 */
template<class Context, class LocalPositionType>
class LocalPointSourceContext : public Context
{
public:
    using LocalPosition = LocalPositionType;

    LocalPointSourceContext() = default;
    LocalPointSourceContext(const Context& ctx, LocalPosition localPosition)
    : Context(ctx)
    , localPosition_(std::move(localPosition))
    {}

    const LocalPosition& localPosition() const
    { return localPosition_; }

private:
    LocalPosition localPosition_;
};

/*!
 * \ingroup Core
 * \brief Combines point source maps with a user-provided evaluation function.
 *
 * A point source acts on a local dof with a control volume if it lies in its sub-control volume (`contexts(elemDisc, scv)`),
 * and on a local dof without control volume, which includes all local dofs of finite element discretizations,
 * with the value of its basis function at the point (`contexts(elemDisc, element)`).
 *
 * \tparam Context the point source context type (e.g. PointSourceContext)
 * \tparam Function callable with signature
 *         `Values(const ElementDiscretization&, const ElementVariables&, const ScvPointSourceContext<Context>&)`
 *         and with `const LocalPointSourceContext<Context, LocalPosition>&` as last argument if the discretization has
 *         local dofs without control volume
 * \tparam LocalPosition the position type in local element coordinates
 */
template<class Context, class Function, class LocalPosition = typename Context::GlobalPosition>
class PointSources
{
    using ScvContext = ScvPointSourceContext<Context>;
    using ElementContext = LocalPointSourceContext<Context, LocalPosition>;
    using Scalar = typename Context::GlobalPosition::value_type;

public:
    //! The contexts of the sub-control volumes, keyed by element index and index of the scv in the element
    using MapType = std::map<std::pair<std::size_t, std::size_t>, std::vector<ScvContext>>;
    //! The contexts of the elements, keyed by element index
    using ElementMapType = std::map<std::size_t, std::vector<ElementContext>>;

    PointSources() = default;
    PointSources(MapType map, Function f)
    : map_(std::move(map)), f_(std::move(f))
    {}

    PointSources(MapType map, ElementMapType elementMap, Function f)
    : map_(std::move(map)), elementMap_(std::move(elementMap)), f_(std::move(f))
    {}

    //! Returns whether no point sources have been registered.
    bool empty() const { return map_.empty() && elementMap_.empty(); }

    //! Return the context range of a sub-control volume.
    template<class ElementDiscretization>
    std::span<const ScvContext> contexts(const ElementDiscretization& elemDisc,
                                         const typename ElementDiscretization::SubControlVolume& scv) const
    {
        const auto eIdx = elemDisc.gridDiscretization().elementMapper().index(elemDisc.element());
        const auto scvIdx = scv.indexInElement();
        const auto it = map_.find({eIdx, scvIdx});
        return (it != map_.end()) ? std::span<const ScvContext>(it->second)
                                  : std::span<const ScvContext>{};
    }

    /*!
     * \brief Return the context range of the whole element, used for the local dofs without control volume.
     * \note Each context is weighted such that the point is counted once in total, also if it lies
     *       in several elements. This holds for the contexts of the sub-control volumes as well.
     * \note The range is empty for discretizations whose local dofs all have a control volume.
     */
    template<class ElementDiscretization>
    std::span<const ElementContext> contexts(const ElementDiscretization& elemDisc,
                                             const typename ElementDiscretization::Element& element) const
    {
        const auto eIdx = elemDisc.gridDiscretization().elementMapper().index(element);
        const auto it = elementMap_.find(eIdx);
        return (it != elementMap_.end()) ? std::span<const ElementContext>(it->second)
                                         : std::span<const ElementContext>{};
    }

    //! Evaluate the source for a given context, divided by the number of embeddings of the context.
    template<class ElementDiscretization, class ElementVariables, class Ctx>
    auto eval(const ElementDiscretization& elemDisc,
              const ElementVariables& elemVars,
              const Ctx& context) const
    {
        auto values = f_(elemDisc, elemVars, context);
        values /= static_cast<Scalar>(context.embeddings());
        return values;
    }

private:
    MapType map_;
    ElementMapType elementMap_;
    Function f_;
};

#ifndef DOXYGEN
namespace PointSourceImpl {

template<class ScvMap, class ElementMap>
struct PointSourceMaps
{
    ScvMap scvMap;
    ElementMap elementMap;
};

template<class GridDiscretization, class Context>
auto makePointSourceMaps(const GridDiscretization& gridDiscretization,
                         const std::vector<Context>& contexts)
{
    using Element = typename GridDiscretization::GridView::template Codim<0>::Entity;
    using LocalPosition = typename Element::Geometry::LocalCoordinate;
    using ScvContext = ScvPointSourceContext<Context>;
    using ElementContext = LocalPointSourceContext<Context, LocalPosition>;
    using ElemDisc = typename GridDiscretization::LocalView;

    PointSourceMaps<std::map<std::pair<std::size_t, std::size_t>, std::vector<ScvContext>>,
                    std::map<std::size_t, std::vector<ElementContext>>> maps;

    constexpr bool locatedContexts = requires (const Context& c) { c.elementIndex(); };
    constexpr bool hasScvs = Concepts::FVElementDiscretization<ElemDisc>
                             && DiscretizationMethods::isCVFE<typename GridDiscretization::DiscretizationMethod>;
    using ScvIndices = Dune::ReservedVector<std::size_t, Dumux::Detail::LocalDofs::maxNumLocalDofs<ElemDisc>()>;

    auto elemDisc = localView(gridDiscretization);
    for (const auto& ctx : contexts)
    {
        // the elements that contain the point, for schemes with control volumes with the sub-control volumes containing it
        std::vector<std::pair<std::size_t, ScvIndices>> located;
        if constexpr (locatedContexts)
            located.emplace_back(ctx.elementIndex(), ScvIndices{});
        else
            for (const auto eIdx : intersectingEntities(ctx.position(), gridDiscretization.boundingBoxTree()))
                located.emplace_back(eIdx, ScvIndices{});

        if constexpr (hasScvs)
        {
            for (auto& [eIdx, scvIndices] : located)
            {
                elemDisc.bindElement(gridDiscretization.element(eIdx));
                scvIndices = scvsContainingPoint<Dumux::Detail::LocalDofs::maxNumLocalDofs<ElemDisc>()>(elemDisc, ctx.position());
                if constexpr (locatedContexts)
                    if (scvIndices.empty())
                        DUNE_THROW(Dune::InvalidStateException,
                                   "Point source at " << ctx.position() << " is not inside any sub control volume of element " << eIdx);
            }
            std::erase_if(located, [](const auto& l){ return l.second.empty(); });
        }

        if (located.empty())
            continue;

        // the number of contexts the point is shared by
        const std::size_t elementEmbeddings = locatedContexts ? ctx.embeddings() : ctx.embeddings()*located.size();
        for (const auto& [eIdx, scvIndices] : located)
        {
            if constexpr (Dumux::Detail::LocalDofs::hasNonCVLocalDofsInterface<ElemDisc>())
            {
                const auto element = gridDiscretization.element(eIdx);
                auto elementEntry = ElementContext(ctx, element.geometry().local(ctx.position()));
                elementEntry.setEmbeddings(elementEmbeddings);
                maps.elementMap[eIdx].push_back(std::move(elementEntry));
            }

            if constexpr (hasScvs)
            {
                for (const auto scvIdx : scvIndices)
                {
                    auto entry = ctx;
                    entry.setEmbeddings(elementEmbeddings * scvIndices.size());
                    maps.scvMap[{eIdx, scvIdx}].emplace_back(entry, scvIdx);
                }
            }
            else if constexpr (GridDiscretization::discMethod == DiscretizationMethods::cctpfa
                            || GridDiscretization::discMethod == DiscretizationMethods::ccmpfa)
            {
                auto entry = ctx;
                entry.setEmbeddings(elementEmbeddings);
                maps.scvMap[{eIdx, std::size_t(0)}].emplace_back(entry, std::size_t(0));
            }
            else if constexpr (!Concepts::FEElementDiscretization<ElemDisc>)
                DUNE_THROW(Dune::NotImplemented,
                           "makePointSourceMap for discretization method " << GridDiscretization::discMethod);
        }
    }
    return maps;
}

} // end namespace PointSourceImpl
#endif

/*!
 * \ingroup Core
 * \brief Build the point source map of the sub-control volumes from a list of point source contexts.
 *
 * Returns `std::map<std::pair<elementIndex, scvIndex>, std::vector<ScvPointSourceContext<Context>>>`.
 * The scv index and embeddings count are set automatically based on the
 * bounding-box-tree intersection of each context's position with the grid.
 * The map is empty for finite element discretizations, which have no sub-control volumes.
 * If the context provides `elementIndex()` (see ElementIndexPointSourceContext), the context is
 * assigned to that element without searching the grid, and its embeddings are the number of
 * elements sharing the point.
 *
 * \tparam GridDiscretization the grid discretization type
 * \tparam Context the point source context type; must provide `position()` and
 *         `setEmbeddings()`
 */
template<class GridDiscretization, class Context>
auto makePointSourceMap(const GridDiscretization& gridDiscretization,
                        std::initializer_list<Context> contexts)
{ return makePointSourceMap(gridDiscretization, std::vector<Context>(contexts)); }

template<class GridDiscretization, class Context>
auto makePointSourceMap(const GridDiscretization& gridDiscretization,
                        const std::vector<Context>& contexts)
{ return PointSourceImpl::makePointSourceMaps(gridDiscretization, contexts).scvMap; }

/*!
 * \ingroup Core
 * \brief Build the point source map of the elements from a list of point source contexts.
 *
 * Returns `std::map<elementIndex, std::vector<LocalPointSourceContext<Context, LocalPosition>>>`.
 * The position in the local coordinates of the element is computed once. The embeddings count is
 * the number of elements that contain the point, such that the weights of its contexts sum to one.
 * The location of the point is determined as in makePointSourceMap.
 * The map is empty for discretizations whose local dofs all have a control volume.
 *
 * \tparam GridDiscretization the grid discretization type
 * \tparam Context the point source context type; must provide `position()` and
 *         `setEmbeddings()`
 */
template<class GridDiscretization, class Context>
auto makeElementPointSourceMap(const GridDiscretization& gridDiscretization,
                               std::initializer_list<Context> contexts)
{ return makeElementPointSourceMap(gridDiscretization, std::vector<Context>(contexts)); }

template<class GridDiscretization, class Context>
auto makeElementPointSourceMap(const GridDiscretization& gridDiscretization,
                               const std::vector<Context>& contexts)
{ return PointSourceImpl::makePointSourceMaps(gridDiscretization, contexts).elementMap; }

/*!
 * \ingroup Core
 * \brief Build a PointSources object from a list of contexts and an evaluation function.
 *
 * \tparam GridDiscretization the grid discretization type
 * \tparam Context the point source context type
 * \tparam Function callable with signature
 *         `Values(const ElementDiscretization&, const ElementVariables&, const ScvPointSourceContext<Context>&)`
 *         and with `const LocalPointSourceContext<Context, LocalPosition>&` as last argument if the discretization has
 *         local dofs without control volume
 */
template<class GridDiscretization, class Context, class Function>
auto makePointSources(const GridDiscretization& gridDiscretization,
                      std::initializer_list<Context> contexts,
                      Function f)
{ return makePointSources(gridDiscretization, std::vector<Context>(contexts), std::move(f)); }

template<class GridDiscretization, class Context, class Function>
auto makePointSources(const GridDiscretization& gridDiscretization,
                      const std::vector<Context>& contexts,
                      Function f)
{
    using Element = typename GridDiscretization::GridView::template Codim<0>::Entity;
    using LocalPosition = typename Element::Geometry::LocalCoordinate;
    auto maps = PointSourceImpl::makePointSourceMaps(gridDiscretization, contexts);
    return PointSources<Context, Function, LocalPosition>(
        std::move(maps.scvMap), std::move(maps.elementMap), std::move(f)
    );
}

} // end namespace Dumux::Experimental

#endif
