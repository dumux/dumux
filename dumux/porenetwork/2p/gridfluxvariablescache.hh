// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PNMTwoPModel
 * \brief Global flux variable cache
 */
#ifndef DUMUX_PNM_2P_GRID_FLUXVARSCACHE_HH
#define DUMUX_PNM_2P_GRID_FLUXVARSCACHE_HH

#include <type_traits>
#include <dune/common/std/type_traits.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/cvfe/gridfluxvariablescache.hh>
#include "elementfluxvariablescache.hh"
#include "invasionstate.hh"

namespace Dumux::PoreNetwork {

namespace Detail {

template<class P>
using ProblemInvasionState = decltype(std::declval<const P&>().invasionState());

//! whether the problem provides the invasion state, i.e. is its owner
template<class P>
constexpr inline bool problemProvidesInvasionState()
{ return Dune::Std::is_detected<ProblemInvasionState, P>::value; }

//! the invasion state as provided by the problem, a reference to the state owned by it
template<class P>
struct InvasionStateFromProblem
{ using type = ProblemInvasionState<P>; };

//! the invasion state as selected by the traits, owned by the grid flux variables cache
template<class Traits>
struct InvasionStateFromTraits
{ using type = typename Traits::InvasionState; };

/*!
 * \brief The type of the invasion state held by the grid flux variables cache.
 * \note Only one of the two is ever named: a reference to the problem's state if it provides one,
 *       the type selected by the traits otherwise, which is the deprecated way of choosing it.
 */
template<class Problem, class Traits>
using InvasionStateType = typename std::conditional_t<problemProvidesInvasionState<Problem>(),
                                                      InvasionStateFromProblem<Problem>,
                                                      InvasionStateFromTraits<Traits>>::type;

} // end namespace Detail

/*!
 * \ingroup PNMTwoPModel
 * \brief Flux variable caches traits
 *
 * \tparam P The problem
 * \tparam FVC The flux variables cache
 * \tparam IS The invasion state.
 *             \deprecated Selecting the invasion state here is deprecated, let the problem own it
 *             instead by providing invasionState(). This parameter and the InvasionState member it
 *             defines will be removed after 3.12.
 */
template<class P, class FVC, class IS = TwoPInvasionState<P>>
struct PNMTwoPDefaultGridFVCTraits
{
    using Problem = P;
    using FluxVariablesCache = FVC;

    //! \deprecated Let the problem own the invasion state instead. Will be removed after 3.12.
    using InvasionState [[deprecated("Let the problem own the invasion state instead. Will be removed after 3.12.")]] = IS;

    template<class GridFluxVariablesCache, bool cachingEnabled>
    using LocalView = PNMTwoPElementFluxVariablesCache<GridFluxVariablesCache, cachingEnabled>;
};

/*!
 * \ingroup PNMTwoPModel
 * \brief Flux variable caches on a gridview
 * \note The class is specialized for a version with and without grid caching
 */
template<class Problem,
         class FluxVariablesCache,
         bool cachingEnabled,
         class Traits>
class PNMTwoPGridFluxVariablesCache;

/*!
 * \ingroup PNMTwoPModel
 * \brief The grid flux variables cache for the two-phase PNM hodling the invasion state of the throats
 * \note The flux caches of the gridview are stored which is memory intensive but faster
 */
template<class P, class FVC, class Traits>
class PNMTwoPGridFluxVariablesCache<P, FVC, true, Traits>
{
    using Problem = typename Traits::Problem;
    using ThisType = PNMTwoPGridFluxVariablesCache<P, FVC, true, Traits>;
    using InvasionStateType = Detail::InvasionStateType<Problem, Traits>;

public:
    //! export the flux variable cache type
    using FluxVariablesCache = typename Traits::FluxVariablesCache;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = true;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    PNMTwoPGridFluxVariablesCache(const Problem& problem)
    : problemPtr_(&problem)
    , invasionState_(makeInvasionState_(problem)) {}

    template<class GridGeometry, class GridVolumeVariables, class SolutionVector>
    void update(const GridGeometry& gridGeometry,
                const GridVolumeVariables& gridVolVars,
                const SolutionVector& sol,
                bool forceUpdate = true)
    {
        fluxVarsCache_.resize(gridGeometry.gridView().size(0));
        auto fvGeometry = localView(gridGeometry);
        auto elemVolVars = localView(gridVolVars);
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            auto eIdx = gridGeometry.elementMapper().index(element);

            // bind the geometries and volume variables to the element (all the elements in stencil)
            fvGeometry.bind(element);
            elemVolVars.bind(element, fvGeometry, sol);

            for (auto&& scvf : scvfs(fvGeometry))
                cache(eIdx, scvf.index()).update(problem(), element, fvGeometry, elemVolVars, scvf, invasionState().invaded(element));
        }
    }

    template<class FVElementGeometry, class ElementVolumeVariables>
    void updateElement(const typename FVElementGeometry::Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars)
    {
        if constexpr (FluxVariablesCache::isSolDependent)
        {
            const auto eIdx = fvGeometry.gridGeometry().elementMapper().index(element);
            for (const auto& scvf : scvfs(fvGeometry))
                cache(eIdx, scvf.index()).update(problem(), element, fvGeometry, elemVolVars, scvf, invasionState().invaded(element));
        }
    }

    const Problem& problem() const
    { return *problemPtr_; }

    // access operator
    const FluxVariablesCache& cache(std::size_t eIdx, [[maybe_unused]] std::size_t scvfIdx) const
    { return fluxVarsCache_[eIdx]; }

    // access operator
    FluxVariablesCache& cache(std::size_t eIdx, [[maybe_unused]] std::size_t scvfIdx)
    { return fluxVarsCache_[eIdx]; }

    /*!
     * \brief Returns the invasion state of the throats
     * \note If the problem provides an invasion state, it is the owner and this only refers to it.
     *       Otherwise the state is owned here, as it always was.
     */
    decltype(auto) invasionState()
    { return (invasionState_); }

    //! \copydoc invasionState()
    decltype(auto) invasionState() const
    { return (invasionState_); }

private:
    //! Refer to the problem's invasion state if it provides one, otherwise construct our own
    static InvasionStateType makeInvasionState_(const Problem& problem)
    {
        if constexpr (Detail::problemProvidesInvasionState<Problem>())
            return problem.invasionState();
        else
            return InvasionStateType(problem);
    }

    const Problem* problemPtr_;
    std::vector<FluxVariablesCache> fluxVarsCache_;
    InvasionStateType invasionState_;
};

/*!
 * \ingroup PNMTwoPModel
 * \brief The grid flux variables cache for the two-phase PNM hodling the invasion state of the throats
 * \note The flux caches of the gridview are stored which is memory intensive but faster
 */
template<class P, class FVC, class Traits>
class PNMTwoPGridFluxVariablesCache<P, FVC, false, Traits>
{
    using Problem = typename Traits::Problem;
    using ThisType = PNMTwoPGridFluxVariablesCache<P, FVC, false, Traits>;
    using InvasionStateType = Detail::InvasionStateType<Problem, Traits>;

    public:
    //! export the flux variable cache type
    using FluxVariablesCache = typename Traits::FluxVariablesCache;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = false;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    PNMTwoPGridFluxVariablesCache(const Problem& problem)
    : problemPtr_(&problem)
    , invasionState_(makeInvasionState_(problem)) {}

    template<class GridGeometry, class GridVolumeVariables, class SolutionVector>
    void update(const GridGeometry& gridGeometry,
                const GridVolumeVariables& gridVolVars,
                const SolutionVector& sol,
                bool forceUpdate = true) {}

    const Problem& problem() const
    { return *problemPtr_; }

    /*!
     * \brief Returns the invasion state of the throats
     * \note If the problem provides an invasion state, it is the owner and this only refers to it.
     *       Otherwise the state is owned here, as it always was.
     */
    decltype(auto) invasionState()
    { return (invasionState_); }

    //! \copydoc invasionState()
    decltype(auto) invasionState() const
    { return (invasionState_); }

private:
    //! Refer to the problem's invasion state if it provides one, otherwise construct our own
    static InvasionStateType makeInvasionState_(const Problem& problem)
    {
        if constexpr (Detail::problemProvidesInvasionState<Problem>())
            return problem.invasionState();
        else
            return InvasionStateType(problem);
    }

    const Problem* problemPtr_;
    InvasionStateType invasionState_;
};
} // end namespace Dumux::PoreNetwork

#endif
