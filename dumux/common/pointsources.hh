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
#include <dumux/geometry/intersectspointgeometry.hh>
#include <dumux/geometry/intersectingentities.hh>
#include <dumux/discretization/method.hh>

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
 * \brief Combines a point source map with a user-provided evaluation function.
 *
 * \tparam Context the point source context type (e.g. PointSourceContext)
 * \tparam Function callable with signature
 *         `Values(const ElementDiscretization&, const ElementVariables&, const ScvPointSourceContext<Context>&)`
 */
template<class Context, class Function>
class PointSources
{
    using ScvContext = ScvPointSourceContext<Context>;
    using Scalar = typename Context::GlobalPosition::value_type;

public:
    using MapType = std::map<std::pair<std::size_t, std::size_t>, std::vector<ScvContext>>;

    PointSources() = default;
    PointSources(MapType map, Function f)
    : map_(std::move(map)), f_(std::move(f))
    {}

    //! Returns whether no point sources have been registered.
    bool empty() const { return map_.empty(); }

    //! Return the context range of a sub-control volume.
    template<class ElementDiscretization>
    std::span<const ScvContext> contexts(const ElementDiscretization& elemDisc,
                                         const typename ElementDiscretization::SubControlVolume& scv) const
    {
        const auto eIdx = elemDisc.gridGeometry().elementMapper().index(elemDisc.element());
        const auto scvIdx = scv.indexInElement();
        const auto it = map_.find({eIdx, scvIdx});
        return (it != map_.end()) ? std::span<const ScvContext>(it->second)
                                  : std::span<const ScvContext>{};
    }

    //! Return the context range of the whole element.
    template<class ElementDiscretization>
    auto contexts(const ElementDiscretization& elemDisc,
                  const typename ElementDiscretization::Element& element) const
    {
        return scvs(elemDisc)
             | std::views::transform([this, &elemDisc](const auto& scv)
               { return this->contexts(elemDisc, scv); })
             | std::views::join;
    }

    //! Evaluate the source for a given context.
    template<class ElementDiscretization, class ElementVariables>
    auto eval(const ElementDiscretization& elemDisc,
              const ElementVariables& elemVars,
              const ScvContext& context) const
    {
        auto values = f_(elemDisc, elemVars, context);
        values /= static_cast<Scalar>(context.embeddings());
        return values;
    }

private:
    MapType map_;
    Function f_;
};

/*!
 * \ingroup Core
 * \brief Build a point source map from a list of point source contexts.
 *
 * Returns `std::map<std::pair<elementIndex, scvIndex>, std::vector<ScvPointSourceContext<Context>>>`.
 * The scv index and embeddings count are set automatically based on the
 * bounding-box-tree intersection of each context's position with the grid.
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
{
    using ScvContext = ScvPointSourceContext<Context>;
    std::map<std::pair<std::size_t, std::size_t>, std::vector<ScvContext>> map;
    const auto& boundingBoxTree = gridDiscretization.boundingBoxTree();

    for (const auto& ctx : contexts)
    {
        const auto entities = intersectingEntities(ctx.position(), boundingBoxTree);
        if (entities.empty())
            continue;

        static constexpr bool isCVFE = DiscretizationMethods::isCVFE<typename GridDiscretization::DiscretizationMethod>;
        if constexpr (isCVFE)
        {
            auto elemDisc = localView(gridDiscretization);
            for (const auto eIdx : entities)
            {
                const auto element = boundingBoxTree.entitySet().entity(eIdx);
                elemDisc.bindElement(element);

                using ElemDisc = typename GridDiscretization::LocalView;
                Dune::ReservedVector<std::size_t, Dumux::Detail::LocalDofs::maxNumLocalDofs<ElemDisc>()> scvIndices;
                for (const auto& scv : scvs(elemDisc))
                    if (intersectsPointGeometry(ctx.position(), elemDisc.geometry(scv)))
                        scvIndices.push_back(scv.indexInElement());

                for (const auto scvIdx : scvIndices)
                {
                    auto entry = ctx;
                    entry.setEmbeddings(ctx.embeddings() * entities.size() * scvIndices.size());
                    map[{eIdx, scvIdx}].emplace_back(entry, scvIdx);
                }
            }
        }
        else if constexpr (GridDiscretization::discMethod == DiscretizationMethods::cctpfa
                        || GridDiscretization::discMethod == DiscretizationMethods::ccmpfa)
        {
            for (const auto eIdx : entities)
            {
                auto entry = ctx;
                entry.setEmbeddings(ctx.embeddings() * entities.size());
                map[{eIdx, std::size_t(0)}].emplace_back(entry, std::size_t(0));
            }
        }
        else
            DUNE_THROW(Dune::NotImplemented,
                       "makePointSourceMap for discretization method " << GridDiscretization::discMethod);
    }
    return map;
}

/*!
 * \ingroup Core
 * \brief Build a PointSources object from a list of contexts and an evaluation function.
 *
 * \tparam GridDiscretization the grid discretization type
 * \tparam Context the point source context type
 * \tparam Function callable with signature
 *         `Values(const ElementDiscretization&, const ElementVariables&, const ScvPointSourceContext<Context>&)`
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
    return PointSources<Context, Function>(
        makePointSourceMap(gridDiscretization, contexts),
        std::move(f)
    );
}

} // end namespace Dumux::Experimental

#endif
