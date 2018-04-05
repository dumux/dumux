// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup StaggeredDiscretization
 * \copydoc Dumux::StaggeredGridFluxVariablesCache
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_GRID_FLUXVARSCACHE_HH
#define DUMUX_DISCRETIZATION_STAGGERED_GRID_FLUXVARSCACHE_HH

//! make the local view function available whenever we use this class
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/staggered/elementfluxvariablescache.hh>

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \brief Traits class to be used for the StaggeredGridVFluxVariablesCache.
 * \tparam P The problem type
 * \tparam FVC The flux variables cache type
 */
template<class P, class FVC>
struct StaggeredDefaultGridFluxVariablesCacheTraits
{
    using Problem = P;
    using FluxVariablesCache = FVC;

    template<class GridFluxVariablesCache, bool cachingEnabled>
    using LocalView = StaggeredElementFluxVariablesCache<GridFluxVariablesCache, cachingEnabled>;
};

/*!
 * \ingroup StaggeredDiscretization
 * \brief Flux variables cache class for staggered models
 */
template<class Problem,
         class FluxVariablesCache,
         bool cachingEnabled = false,
         class Traits = StaggeredDefaultGridFluxVariablesCacheTraits<Problem, FluxVariablesCache> >
class StaggeredGridFluxVariablesCache;

/*!
 * \ingroup StaggeredDiscretization
 * \brief Flux variables cache class for staggered models.
          Specialization in case of storing the flux cache.
 */
template<class P, class FVC, class Traits>
class StaggeredGridFluxVariablesCache<P, FVC, true, Traits>
{
    using ThisType = StaggeredGridFluxVariablesCache<P, FVC, true, Traits>;
    using Problem = typename Traits::Problem;

public:
    //! export the flux variable cache type
    using FluxVariablesCache = typename Traits::FluxVariablesCache;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = true;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    StaggeredGridFluxVariablesCache(const Problem& problem) : problemPtr_(&problem) {}

    // When global caching is enabled, precompute transmissibilities and stencils for all the scv faces
    template<class FVGridGeometry, class GridVolumeVariables, class SolutionVector>
    void update(const FVGridGeometry& fvGridGeometry,
                const GridVolumeVariables& gridVolVars,
                const SolutionVector& sol,
                bool forceUpdate = false)
    {
        // fluxVarsCache_.resize(fvGridGeometry.numScvf());
        // for (const auto& element : elements(fvGridGeometry.gridView()))
        // {
        //     // Prepare the geometries within the elements of the stencil
        //     auto fvGeometry = localView(fvGridGeometry);
        //     fvGeometry.bind(element);
        //
        //     auto elemVolVars = localView(gridVolVars);
        //     elemVolVars.bind(element, fvGeometry, sol);
        //
        //     for (auto&& scvf : scvfs(fvGeometry))
        //     {
        //         fluxVarsCache_[scvf.index()].update(problem, element, fvGeometry, elemVolVars, scvf);
        //     }
        // }
    }

    const Problem& problem() const
    { return *problemPtr_; }

    // access operators in the case of caching
    const FluxVariablesCache& operator [](std::size_t scvfIdx) const
    { return fluxVarsCache_[scvfIdx]; }

    FluxVariablesCache& operator [](std::size_t scvfIdx)
    { return fluxVarsCache_[scvfIdx]; }

private:
    const Problem* problemPtr_;

    std::vector<FluxVariablesCache> fluxVarsCache_;
    std::vector<std::size_t> globalScvfIndices_;
};

/*!
 * \ingroup StaggeredDiscretization
 * \brief Flux variables cache class for staggered models.
          Specialization in case of not storing the flux cache.
 */
template<class P, class FVC, class Traits>
class StaggeredGridFluxVariablesCache<P, FVC, false, Traits>
{
    using ThisType = StaggeredGridFluxVariablesCache<P, FVC, false, Traits>;
    using Problem = typename Traits::Problem;

public:
    //! export the flux variable cache type
    using FluxVariablesCache = typename Traits::FluxVariablesCache;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = true;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    // When global flux variables caching is disabled, we don't need to update the cache
    void update(Problem& problem)
    { problemPtr_ = &problem; }

private:

    const Problem& problem_() const
    { return *problemPtr_; }

    const Problem* problemPtr_;
};

} // end namespace Dumux

#endif
