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
 * \ingroup StaggeredDiscretization
 * \copydoc Dumux::StaggeredGridFluxVariablesCache
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_GRID_FLUXVARSCACHE_HH
#define DUMUX_DISCRETIZATION_STAGGERED_GRID_FLUXVARSCACHE_HH

// make the local view function available whenever we use this class
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/staggered/elementfluxvariablescache.hh>

#include <dumux/freeflow/staggeredupwindmethods.hh>

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \brief Traits class to be used for the StaggeredGridVFluxVariablesCache.
 * \tparam P The problem type
 * \tparam FVC The flux variables cache type
 */
template<class P, class FVC, class FVCF, int upwOrder>
struct StaggeredDefaultGridFluxVariablesCacheTraits
{
    using Problem = P;
    using FluxVariablesCache = FVC;
    using FluxVariablesCacheFiller = FVCF;

    template<class GridFluxVariablesCache, bool cachingEnabled>
    using LocalView = StaggeredElementFluxVariablesCache<GridFluxVariablesCache, cachingEnabled>;
    static constexpr int upwindSchemeOrder = upwOrder;
};

/*!
 * \ingroup StaggeredDiscretization
 * \brief Flux variables cache class for staggered models
 */
template<class Problem,
         class FluxVariablesCache,
         class FluxVariablesCacheFiller,
         bool EnableGridFluxVariablesCache = false,
         int upwindSchemeOrder = 1,
         class Traits = StaggeredDefaultGridFluxVariablesCacheTraits<Problem, FluxVariablesCache, FluxVariablesCacheFiller, upwindSchemeOrder>>
class StaggeredGridFluxVariablesCache;

/*!
 * \ingroup StaggeredDiscretization
 * \brief Flux variables cache class for staggered models.
          Specialization in case of storing the flux cache.
 */
template<class P, class FVC, class FVCF, int upwindSchemeOrder, class TheTraits>
class StaggeredGridFluxVariablesCache<P, FVC, FVCF, true, upwindSchemeOrder, TheTraits>
{
    using Problem = typename TheTraits::Problem;
    using ThisType = StaggeredGridFluxVariablesCache<P, FVC, FVCF, true, upwindSchemeOrder, TheTraits>;

    //!  the flux variable cache filler type
    using FluxVariablesCacheFiller = typename TheTraits::FluxVariablesCacheFiller;
public:
    //! the flux var cache traits
    using Traits = TheTraits;

    //! export the flux variable cache type
    using FluxVariablesCache = typename Traits::FluxVariablesCache;
    using Scalar = typename FluxVariablesCache::Scalar;

    static constexpr bool useHigherOrder = upwindSchemeOrder > 1;

    //! export upwind scheme
    using UpwindScheme = StaggeredUpwindMethods<Scalar, upwindSchemeOrder>;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = true;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    StaggeredGridFluxVariablesCache(const Problem& problem)
    : problemPtr_(&problem)
    , staggeredUpwindMethods_(problem.paramGroup())
    {}

    // When global caching is enabled, precompute transmissibilities and stencils for all the scv faces
    template<class GridGeometry, class GridVolumeVariables, class SolutionVector>
    void update(const GridGeometry& gridGeometry,
                const GridVolumeVariables& gridVolVars,
                const SolutionVector& sol,
                bool forceUpdate = false)
    {
        // only do the update if fluxes are solution dependent or if update is forced
        if (FluxVariablesCacheFiller::isSolDependent || forceUpdate)
        {
            // instantiate helper class to fill the caches
            // FluxVariablesCacheFiller filler(problem()); TODO: use proper ctor
            FluxVariablesCacheFiller filler(problem());

            fluxVarsCache_.resize(gridGeometry.numScvf());
            for (const auto& element : elements(gridGeometry.gridView()))
            {
                // Prepare the geometries within the elements of the stencil
                auto fvGeometry = localView(gridGeometry);
                fvGeometry.bind(element);

                auto elemVolVars = localView(gridVolVars);
                elemVolVars.bind(element, fvGeometry, sol);

                for (auto&& scvf : scvfs(fvGeometry))
                {
                    filler.fill(*this, fluxVarsCache_[scvf.index()], element, fvGeometry, elemVolVars, scvf, forceUpdate);
                }
            }
        }
    }

    //! Return the StaggeredUpwindMethods
    const UpwindScheme& staggeredUpwindMethods() const
    {
        return staggeredUpwindMethods_;
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
    UpwindScheme staggeredUpwindMethods_;

    std::vector<FluxVariablesCache> fluxVarsCache_;
    std::vector<std::size_t> globalScvfIndices_;
};

/*!
 * \ingroup StaggeredDiscretization
 * \brief Flux variables cache class for staggered models.
          Specialization in case of not storing the flux cache.
 */
template<class P, class FVC, class FVCF, int upwindSchemeOrder, class TheTraits>
class StaggeredGridFluxVariablesCache<P, FVC, FVCF, false, upwindSchemeOrder, TheTraits>
{
    using Problem = typename TheTraits::Problem;
    using ThisType = StaggeredGridFluxVariablesCache<P, FVC, FVCF, false, upwindSchemeOrder, TheTraits>;

    //!  the flux variable cache filler type
    using FluxVariablesCacheFiller = typename TheTraits::FluxVariablesCacheFiller;
public:
    //! the flux var cache traits
    using Traits = TheTraits;

    //! export the flux variable cache type
    using FluxVariablesCache = typename Traits::FluxVariablesCache;

    //! the scalar type
    using Scalar = typename FluxVariablesCache::Scalar;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = false;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    //! export upwind scheme
    using UpwindScheme = StaggeredUpwindMethods<Scalar, upwindSchemeOrder>;

    StaggeredGridFluxVariablesCache(const Problem& problem)
    : problemPtr_(&problem)
    , staggeredUpwindMethods_(problem.paramGroup())
      {}

    // When global caching is enabled, precompute transmissibilities and stencils for all the scv faces
    template<class GridGeometry, class GridVolumeVariables, class SolutionVector>
    void update(const GridGeometry& gridGeometry,
                const GridVolumeVariables& gridVolVars,
                const SolutionVector& sol,
                bool forceUpdate = false) {}

    const Problem& problem() const
    { return *problemPtr_; }

    //! Return the UpwindingMethods
    const UpwindScheme& staggeredUpwindMethods() const
    {
        return staggeredUpwindMethods_;
    }

private:
    const Problem* problemPtr_;
    UpwindScheme staggeredUpwindMethods_;
};

} // end namespace Dumux

#endif
