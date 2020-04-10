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
 * \ingroup BoxDiscretization
 * \brief Global flux variable cache
 */
#ifndef DUMUX_DISCRETIZATION_BOX_GRID_FLUXVARSCACHE_HH
#define DUMUX_DISCRETIZATION_BOX_GRID_FLUXVARSCACHE_HH

// make the local view function available whenever we use this class
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/box/elementfluxvariablescache.hh>

namespace Dumux {

/*!
 * \ingroup BoxDiscretization
 * \brief Flux variable caches traits
 */
template<class P, class FVC>
struct BoxDefaultGridFVCTraits
{
    using Problem = P;
    using FluxVariablesCache = FVC;

    template<class GridFluxVariablesCache, bool cachingEnabled>
    using LocalView = BoxElementFluxVariablesCache<GridFluxVariablesCache, cachingEnabled>;
};

/*!
 * \ingroup BoxDiscretization
 * \brief Flux variable caches on a gridview
 * \note The class is specialized for a version with and without grid caching
 */
template<class Problem,
         class FluxVariablesCache,
         bool cachingEnabled = false,
         class Traits = BoxDefaultGridFVCTraits<Problem, FluxVariablesCache> >
class BoxGridFluxVariablesCache;

/*!
 * \ingroup BoxDiscretization
 * \brief Flux variable caches on a gridview with grid caching enabled
 * \note The flux caches of the gridview are stored which is memory intensive but faster
 */
template<class P, class FVC, class Traits>
class BoxGridFluxVariablesCache<P, FVC, true, Traits>
{
    using Problem = typename Traits::Problem;
    using ThisType = BoxGridFluxVariablesCache<P, FVC, true, Traits>;

public:
    //! export the flux variable cache type
    using FluxVariablesCache = typename Traits::FluxVariablesCache;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = true;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    BoxGridFluxVariablesCache(const Problem& problem) : problemPtr_(&problem) {}

    template<class GridGeometry, class GridVolumeVariables, class SolutionVector>
    void update(const GridGeometry& gridGeometry,
                const GridVolumeVariables& gridVolVars,
                const SolutionVector& sol,
                bool forceUpdate = false)
    {
        // Here, we do not do anything unless it is a forced update
        if (forceUpdate)
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
                    cache(eIdx, scvf.index()).update(problem(), element, fvGeometry, elemVolVars, scvf);
            }
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

private:
    // currently bound element
    const Problem* problemPtr_;
    std::vector<std::vector<FluxVariablesCache>> fluxVarsCache_;
};

/*!
 * \ingroup BoxDiscretization
 * \brief Flux variable caches on a gridview with grid caching disabled
 */
template<class P, class FVC, class Traits>
class BoxGridFluxVariablesCache<P, FVC, false, Traits>
{
    using Problem = typename Traits::Problem;
    using ThisType = BoxGridFluxVariablesCache<P, FVC, false, Traits>;

public:
    //! export the flux variable cache type
    using FluxVariablesCache = typename Traits::FluxVariablesCache;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = false;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    BoxGridFluxVariablesCache(const Problem& problem) : problemPtr_(&problem) {}

    template<class GridGeometry, class GridVolumeVariables, class SolutionVector>
    void update(const GridGeometry& gridGeometry,
                const GridVolumeVariables& gridVolVars,
                const SolutionVector& sol,
                bool forceUpdate = false) {}

    const Problem& problem() const
    { return *problemPtr_; }

private:
    const Problem* problemPtr_;
};

} // end namespace Dumux

#endif
