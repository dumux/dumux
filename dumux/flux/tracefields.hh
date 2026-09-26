// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Flux
 * \brief Fields over an element, evaluated lazily from a flux law
 */
#ifndef DUMUX_FLUX_TRACE_FIELDS_HH
#define DUMUX_FLUX_TRACE_FIELDS_HH

#include <memory>
#include <type_traits>
#include <utility>
#include <variant>

#include <dumux/common/concepts/variables_.hh>
#include <dumux/common/typetraits/griddiscretization.hh>
#include <dumux/discretization/localview.hh>

namespace Dumux {

/*!
 * \ingroup Flux
 * \brief The element-local data a flux law needs, gathered in one place.
 *
 * Holds its local views by value and keeps its dependencies alive through shared
 * ownership, so a bound context stays valid however the field it came from was created.
 */
template<class Problem, class GridVariables, class SolutionVector>
class ElementFluxContext
{
    using GridGeometry = std::decay_t<decltype(gridDiscretization(std::declval<const Problem&>()))>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;

    //! Variables are held either as separate volume and flux caches or as one combined cache
    static constexpr bool hasSeparateCaches = Concept::FVGridVariables<GridVariables>;

    static auto makeElemVars(const GridVariables& gridVariables)
    {
        if constexpr (hasSeparateCaches) return localView(gridVariables.curGridVolVars());
        else return localView(gridVariables.curGridVars());
    }

    static auto makeFluxCache(const GridVariables& gridVariables)
    {
        if constexpr (hasSeparateCaches) return localView(gridVariables.gridFluxVarsCache());
        else return std::monostate{};
    }

    using ElemVars = decltype(makeElemVars(std::declval<const GridVariables&>()));
    using FluxCache = decltype(makeFluxCache(std::declval<const GridVariables&>()));

public:
    //! Bind the local views of the discretization and the variables to the given element
    ElementFluxContext(std::shared_ptr<const Problem> problem,
                       std::shared_ptr<const GridVariables> gridVariables,
                       std::shared_ptr<const SolutionVector> x,
                       const Element& element)
    : problem_(std::move(problem))
    , gridVariables_(std::move(gridVariables))
    , x_(std::move(x))
    , element_(element)
    , fvGeometry_(localView(gridDiscretization(*problem_)))
    , elemVars_(makeElemVars(*gridVariables_))
    , fluxCache_(makeFluxCache(*gridVariables_))
    {
        fvGeometry_.bind(element_);
        elemVars_.bind(element_, fvGeometry_, *x_);
        if constexpr (hasSeparateCaches)
            fluxCache_.bind(element_, fvGeometry_, elemVars_);
    }

    //! The problem
    const Problem& problem() const { return *problem_; }
    //! The element the context is bound to
    const Element& element() const { return element_; }
    //! The local view of the discretization, bound to the element
    const auto& fvGeometry() const { return fvGeometry_; }
    //! The local view of the variables, bound to the element
    const auto& elemVolVars() const { return elemVars_; }

    //! The local view of the flux variables cache, the variables themselves if they combine both
    const auto& elemFluxVarsCache() const
    {
        if constexpr (hasSeparateCaches) return fluxCache_;
        else return elemVars_;
    }

private:
    std::shared_ptr<const Problem> problem_;
    std::shared_ptr<const GridVariables> gridVariables_;
    std::shared_ptr<const SolutionVector> x_;
    Element element_;
    typename GridGeometry::LocalView fvGeometry_;
    ElemVars elemVars_;
    FluxCache fluxCache_;
};

/*!
 * \ingroup Flux
 * \brief A field over an element, evaluated lazily from a bound context.
 *
 * The invoker supplies the physics: it receives the bound context and the entity to evaluate
 * at, and calls whatever flux law it wraps. Flux laws differ in the order and even the kind of
 * their arguments, so absorbing that difference here is what lets consumers of this field stay
 * ignorant of the physics. What the field is evaluated at is the invoker's business, a
 * sub-control volume face or an oriented interpolation point alike.
 *
 * \tparam Problem The problem type
 * \tparam GridVariables The grid variables type
 * \tparam SolutionVector The solution vector type
 * \tparam Invoker Callable as invoker(context, entity)
 */
template<class Problem, class GridVariables, class SolutionVector, class Invoker>
class LazyFluxField
{
    using Context = ElementFluxContext<Problem, GridVariables, SolutionVector>;
    using GridGeometry = std::decay_t<decltype(gridDiscretization(std::declval<const Problem&>()))>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;

public:
    class LocalView
    {
    public:
        LocalView(Context context, Invoker invoker)
        : context_(std::move(context)), invoker_(std::move(invoker)) {}

        //! Evaluate the field at the given entity of the bound element
        template<class Entity>
        auto operator()(const Entity& entity) const
        { return invoker_(context_, entity); }

        //! The context the field is evaluated in
        const Context& context() const { return context_; }

    private:
        Context context_;
        Invoker invoker_;
    };

    LazyFluxField(std::shared_ptr<const Problem> problem,
                  std::shared_ptr<const GridVariables> gridVariables,
                  std::shared_ptr<const SolutionVector> x,
                  Invoker invoker)
    : problem_(std::move(problem))
    , gridVariables_(std::move(gridVariables))
    , x_(std::move(x))
    , invoker_(std::move(invoker))
    {}

    //! Bind the field to an element; the result shares ownership of what the field reads
    LocalView bind(const Element& element) const
    { return LocalView{Context{problem_, gridVariables_, x_, element}, invoker_}; }

private:
    std::shared_ptr<const Problem> problem_;
    std::shared_ptr<const GridVariables> gridVariables_;
    std::shared_ptr<const SolutionVector> x_;
    Invoker invoker_;
};

/*!
 * \ingroup Flux
 * \brief A lazily evaluated flux field sharing ownership of what it reads.
 */
template<class Problem, class GridVariables, class SolutionVector, class Invoker>
auto makeLazyFluxField(std::shared_ptr<const Problem> problem,
                       std::shared_ptr<const GridVariables> gridVariables,
                       std::shared_ptr<const SolutionVector> x,
                       Invoker invoker)
{
    return LazyFluxField<Problem, GridVariables, SolutionVector, Invoker>{
        std::move(problem), std::move(gridVariables), std::move(x), std::move(invoker)
    };
}

/*!
 * \ingroup Flux
 * \brief A lazily evaluated flux field that does not own what it reads.
 * \note The caller guarantees that problem, grid variables and solution outlive every object
 *       bound from the field. Prefer makeLazyFluxField, which owns them, wherever the
 *       field may outlive the scope it was created in.
 */
template<class Problem, class GridVariables, class SolutionVector, class Invoker>
auto makeNonOwningFluxField(const Problem& problem,
                            const GridVariables& gridVariables,
                            const SolutionVector& x,
                            Invoker invoker)
{
    const auto alias = [] (const auto& t) {
        using T = std::decay_t<decltype(t)>;
        return std::shared_ptr<const T>(&t, [] (const T*) {});
    };
    return makeLazyFluxField(alias(problem), alias(gridVariables), alias(x), std::move(invoker));
}

} // end namespace Dumux

#endif
