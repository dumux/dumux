// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup PoreNetworkDiscretization
 * \brief Global flux variable cache
 */
#ifndef DUMUX_PNM_2P_GRID_FLUXVARSCACHE_HH
#define DUMUX_PNM_2P_GRID_FLUXVARSCACHE_HH

#include <dumux/common/parameters.hh>
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/box/gridfluxvariablescache.hh>
#include "elementfluxvariablescache.hh"
#include "invasionstate.hh"

namespace Dumux {

/*!
 * \ingroup PoreNetworkDiscretization
 * \brief Flux variable caches traits
 */
template<class P, class FVC, class I, class L>
struct PNMTwoPDefaultGridFVCTraits
{
    using Problem = P;
    using FluxVariablesCache = FVC;
    using Indices = I;
    using Labels = L;

    template<class GridFluxVariablesCache, bool cachingEnabled>
    using LocalView = PNMTwoPElementFluxVariablesCache<GridFluxVariablesCache, cachingEnabled>;
};

/*!
 * \ingroup BoxDiscretization
 * \brief Flux variable caches on a gridview
 * \note The class is specialized for a version with and without grid caching
 */
template<class Problem,
         class FluxVariablesCache,
         bool cachingEnabled,
         class Traits>
class PNMTwoPGridFluxVariablesCache;

/*!
 * \ingroup PoreNetworkDiscretization
 * \brief The grid flux variables cache for the two-phase PNM hodling the invasion state of the throats
 * \note The flux caches of the gridview are stored which is memory intensive but faster
 */
template<class P, class FVC, class Traits>
class PNMTwoPGridFluxVariablesCache<P, FVC, true, Traits>
{
    using Problem = typename Traits::Problem;
    using ThisType = PNMTwoPGridFluxVariablesCache<P, FVC, true, Traits>;
    using InvasionState = PNMTwoPInvasionState<Problem>;

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
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            auto eIdx = gridGeometry.elementMapper().index(element);
            // bind the geometries and volume variables to the element (all the elements in stencil)
            auto fvGeometry = localView(gridGeometry);
            fvGeometry.bind(element);

            auto elemVolVars = localView(gridVolVars);
            elemVolVars.bind(element, fvGeometry, sol);

            fluxVarsCache_[eIdx].resize(fvGeometry.numScvf());
            for (auto&& scvf : scvfs(fvGeometry))
                cache(eIdx, scvf.index()).update(problem(), element, fvGeometry, elemVolVars, scvf, invasionState().invaded(element));
        }
    }

    const Problem& problem() const
    { return *problemPtr_; }

    // access operator
    const FluxVariablesCache& cache(std::size_t eIdx, std::size_t scvfIdx) const
    { return fluxVarsCache_[eIdx][scvfIdx]; }

    // access operator
    FluxVariablesCache& cache(std::size_t eIdx, std::size_t scvfIdx)
    { return fluxVarsCache_[eIdx][scvfIdx]; }

    const InvasionState& invasionState() const
    { return invasionState_; }

    InvasionState& invasionState()
    { return invasionState_; }

    private:
    // currently bound element
    const Problem* problemPtr_;
    std::vector<std::vector<FluxVariablesCache>> fluxVarsCache_;

    InvasionState invasionState_;
};

/*!
 * \ingroup PoreNetworkDiscretization
 * \brief The grid flux variables cache for the two-phase PNM hodling the invasion state of the throats
 * \note The flux caches of the gridview are stored which is memory intensive but faster
 */
template<class P, class FVC, class Traits>
class PNMTwoPGridFluxVariablesCache<P, FVC, false, Traits>
{
    using Problem = typename Traits::Problem;
    using ThisType = PNMTwoPGridFluxVariablesCache<P, FVC, false, Traits>;
    using InvasionState = PNMTwoPInvasionState<Problem>;

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
    // currently bound element
    const Problem* problemPtr_;
    std::vector<std::vector<FluxVariablesCache>> fluxVarsCache_;

    InvasionState invasionState_;
};
} // end namespace Dumux

#endif
