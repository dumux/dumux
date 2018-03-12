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

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \brief Flux variables cache class for staggered models
 */
template<class FVGridGeometry, class Traits, bool EnableGridFluxVariablesCache>
class StaggeredGridFluxVariablesCache;

/*!
 * \ingroup StaggeredDiscretization
 * \brief Flux variables cache class for staggered models.
          Specialization in case of storing the flux cache.
 */
template<class FVGridGeometry, class Traits>
class StaggeredGridFluxVariablesCache<FVGridGeometry, Traits, true>
{
    using ThisType = StaggeredGridFluxVariablesCache<FVGridGeometry, Traits, true>;
    using Problem = typename Traits::Problem;
    using IndexType = typename FVGridGeometry::GridView::IndexSet::IndexType;

public:
    //! export the type of the flux variables cache
    using FluxVariablesCache = typename Traits::FluxVariablesCache;
    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<FVGridGeometry, ThisType, true>;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = true;

    StaggeredGridFluxVariablesCache(const Problem& problem) : problemPtr_(&problem) {}

    // When global caching is enabled, precompute transmissibilities and stencils for all the scv faces
    template<class GridVolumeVariables, class SolutionVector>
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
    const FluxVariablesCache& operator [](IndexType scvfIdx) const
    { return fluxVarsCache_[scvfIdx]; }

    FluxVariablesCache& operator [](IndexType scvfIdx)
    { return fluxVarsCache_[scvfIdx]; }

private:

    const Problem* problemPtr_;

    std::vector<FluxVariablesCache> fluxVarsCache_;
    std::vector<IndexType> globalScvfIndices_;
};

/*!
 * \ingroup StaggeredDiscretization
 * \brief Flux variables cache class for staggered models.
          Specialization in case of not storing the flux cache.
 */
template<class FVGridGeometry, class Traits>
class StaggeredGridFluxVariablesCache<FVGridGeometry, Traits, false>
{
    using ThisType = StaggeredGridFluxVariablesCache<FVGridGeometry, Traits, false>;
    using Problem = typename Traits::Problem;

public:
    //! export the type of the flux variables cache
    using FluxVariablesCache = typename Traits::FluxVariablesCache;
    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<FVGridGeometry, ThisType, false>;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = false;

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
