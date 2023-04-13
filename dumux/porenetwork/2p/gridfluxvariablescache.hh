// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PNMTwoPModel
 * \brief Global flux variable cache
 */
#ifndef DUMUX_PNM_2P_GRID_FLUXVARSCACHE_HH
#define DUMUX_PNM_2P_GRID_FLUXVARSCACHE_HH

#include <dumux/common/parameters.hh>
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/cvfe/gridfluxvariablescache.hh>
#include "elementfluxvariablescache.hh"
#include "invasionstate.hh"

namespace Dumux::PoreNetwork {

/*!
 * \ingroup PNMTwoPModel
 * \brief Flux variable caches traits
 */
template<class P, class FVC, class IS = TwoPInvasionState<P>>
struct PNMTwoPDefaultGridFVCTraits
{
    using Problem = P;
    using FluxVariablesCache = FVC;
    using InvasionState = IS;

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
    using InvasionState = typename Traits::InvasionState;

public:
    //! export the flux variable cache type
    using FluxVariablesCache = typename Traits::FluxVariablesCache;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = true;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    PNMTwoPGridFluxVariablesCache(const Problem& problem)
    : problemPtr_(&problem)
    , invasionState_(problem) {}

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

    const InvasionState& invasionState() const
    { return invasionState_; }

    InvasionState& invasionState()
    { return invasionState_; }

private:
    const Problem* problemPtr_;
    std::vector<FluxVariablesCache> fluxVarsCache_;
    InvasionState invasionState_;
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
    using InvasionState = typename Traits::InvasionState;

    public:
    //! export the flux variable cache type
    using FluxVariablesCache = typename Traits::FluxVariablesCache;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = false;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    PNMTwoPGridFluxVariablesCache(const Problem& problem)
    : problemPtr_(&problem)
    , invasionState_(problem) {}

    template<class GridGeometry, class GridVolumeVariables, class SolutionVector>
    void update(const GridGeometry& gridGeometry,
                const GridVolumeVariables& gridVolVars,
                const SolutionVector& sol,
                bool forceUpdate = true) {}

    const Problem& problem() const
    { return *problemPtr_; }

    const InvasionState& invasionState() const
    { return invasionState_; }

    InvasionState& invasionState()
    { return invasionState_; }

private:
    const Problem* problemPtr_;
    InvasionState invasionState_;
};
} // end namespace Dumux::PoreNetwork

#endif
