// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CVFEDiscretization
 * \brief Global flux variable cache
 */
#ifndef DUMUX_DISCRETIZATION_CVFE_GRID_FLUXVARSCACHE_HH
#define DUMUX_DISCRETIZATION_CVFE_GRID_FLUXVARSCACHE_HH

#include <dumux/parallel/parallel_for.hh>

// make the local view function available whenever we use this class
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/cvfe/elementfluxvariablescache.hh>

namespace Dumux {

/*!
 * \ingroup CVFEDiscretization
 * \brief Flux variable caches traits
 */
template<class P, class FVC>
struct CVFEDefaultGridFVCTraits
{
    using Problem = P;
    using FluxVariablesCache = FVC;

    template<class GridFluxVariablesCache, bool cachingEnabled>
    using LocalView = CVFEElementFluxVariablesCache<GridFluxVariablesCache, cachingEnabled>;
};

/*!
 * \ingroup CVFEDiscretization
 * \brief Flux variable caches on a gridview
 * \note The class is specialized for a version with and without grid caching
 */
template<class Problem,
         class FluxVariablesCache,
         bool cachingEnabled = false,
         class Traits = CVFEDefaultGridFVCTraits<Problem, FluxVariablesCache> >
class CVFEGridFluxVariablesCache;

/*!
 * \ingroup CVFEDiscretization
 * \brief Flux variable caches on a gridview with grid caching enabled
 * \note The flux caches of the gridview are stored which is memory intensive but faster
 */
template<class P, class FVC, class Traits>
class CVFEGridFluxVariablesCache<P, FVC, true, Traits>
{
    using Problem = typename Traits::Problem;
    using ThisType = CVFEGridFluxVariablesCache<P, FVC, true, Traits>;

public:
    //! export the flux variable cache type
    using FluxVariablesCache = typename Traits::FluxVariablesCache;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = true;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    CVFEGridFluxVariablesCache(const Problem& problem) : problemPtr_(&problem) {}

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
            Dumux::parallelFor(gridGeometry.gridView().size(0), [&, &problem = problem()](const std::size_t eIdx)
            {
                // Prepare the geometries within the elements of the stencil
                const auto element = gridGeometry.element(eIdx);
                const auto fvGeometry = localView(gridGeometry).bind(element);
                const auto elemVolVars = localView(gridVolVars).bind(element, fvGeometry, sol);

                // only update shape functions for fluxes if update is forced
                fluxVarsCache_[eIdx].resize(fvGeometry.numScvf());
                for (const auto& scvf : scvfs(fvGeometry))
                    cache(eIdx, scvf.index()).update(problem, element, fvGeometry, elemVolVars, scvf);
            });
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
            fluxVarsCache_[eIdx].resize(fvGeometry.numScvf());
            for (const auto& scvf : scvfs(fvGeometry))
                cache(eIdx, scvf.index()).update(problem(), element, fvGeometry, elemVolVars, scvf);
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
 * \ingroup CVFEDiscretization
 * \brief Flux variable caches on a gridview with grid caching disabled
 */
template<class P, class FVC, class Traits>
class CVFEGridFluxVariablesCache<P, FVC, false, Traits>
{
    using Problem = typename Traits::Problem;
    using ThisType = CVFEGridFluxVariablesCache<P, FVC, false, Traits>;

public:
    //! export the flux variable cache type
    using FluxVariablesCache = typename Traits::FluxVariablesCache;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = false;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    CVFEGridFluxVariablesCache(const Problem& problem) : problemPtr_(&problem) {}

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
