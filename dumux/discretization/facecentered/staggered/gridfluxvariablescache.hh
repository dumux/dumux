// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief Global flux variable cache
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_GRID_FLUXVARSCACHE_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_GRID_FLUXVARSCACHE_HH

#include <dumux/parallel/parallel_for.hh>

// make the local view function available whenever we use this class
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/facecentered/staggered/elementfluxvariablescache.hh>

namespace Dumux {

/*!
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief Flux variable caches traits
 */
template<class P, class FVC, class FVCF>
struct FaceCenteredStaggeredDefaultGridFVCTraits
{
    using Problem = P;
    using FluxVariablesCache = FVC;
    using FluxVariablesCacheFiller = FVCF;

    template<class GridFluxVariablesCache, bool cachingEnabled>
    using LocalView = FaceCenteredStaggeredElementFluxVariablesCache<GridFluxVariablesCache, cachingEnabled>;
};

/*!
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief Flux variable caches on a gridview
 * \note The class is specialized for a version with and without grid caching
 */
template<class Problem,
         class FluxVariablesCache,
         class FluxVariablesCacheFiller,
         bool cachingEnabled = false,
         class Traits = FaceCenteredStaggeredDefaultGridFVCTraits<Problem, FluxVariablesCache, FluxVariablesCacheFiller>>
class FaceCenteredStaggeredGridFluxVariablesCache;

/*!
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief Flux variable caches on a gridview with grid caching enabled
 * \note The flux caches of the gridview are stored which is memory intensive but faster
 */
template<class P, class FVC, class FVCF, class Traits>
class FaceCenteredStaggeredGridFluxVariablesCache<P, FVC, FVCF, true, Traits>
{
    using Problem = typename Traits::Problem;
    using ThisType = FaceCenteredStaggeredGridFluxVariablesCache<P, FVC, FVCF, true, Traits>;
    using FluxVariablesCacheFiller = typename Traits::FluxVariablesCacheFiller;

public:
    //! export the flux variable cache type
    using FluxVariablesCache = typename Traits::FluxVariablesCache;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = true;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    FaceCenteredStaggeredGridFluxVariablesCache(const Problem& problem) : problemPtr_(&problem) {}

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
            FluxVariablesCacheFiller filler(problem());

            fluxVarsCache_.resize(gridGeometry.numScvf());

            Dumux::parallelFor(gridGeometry.gridView().size(0), [&](const std::size_t eIdx)
            {
                // Prepare the geometries within the elements of the stencil
                const auto element = gridGeometry.element(eIdx);
                const auto fvGeometry = localView(gridGeometry).bind(element);
                const auto elemVolVars = localView(gridVolVars).bind(element, fvGeometry, sol);

                for (auto&& scvf : scvfs(fvGeometry))
                {
                    filler.fill(*this, fluxVarsCache_[scvf.index()], element, fvGeometry, elemVolVars, scvf, forceUpdate);
                }
            });
        }
    }

    const Problem& problem() const
    { return *problemPtr_; }

    template<class SubControlVolumeFace>
    const FluxVariablesCache& operator [](const SubControlVolumeFace& scvf) const
    { return fluxVarsCache_[scvf.index()]; }

    template<class SubControlVolumeFace>
    FluxVariablesCache& operator [](const SubControlVolumeFace& scvf)
    { return fluxVarsCache_[scvf.index()]; }

private:
    // currently bound element
    const Problem* problemPtr_;
    std::vector<FluxVariablesCache> fluxVarsCache_;
};

/*!
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief Flux variable caches on a gridview with grid caching disabled
 */
template<class P, class FVC, class FVCF, class Traits>
class FaceCenteredStaggeredGridFluxVariablesCache<P, FVC, FVCF, false, Traits>
{
    using Problem = typename Traits::Problem;
    using ThisType = FaceCenteredStaggeredGridFluxVariablesCache<P, FVC, FVCF, false, Traits>;

public:
    //! export the flux variable cache type
    using FluxVariablesCache = typename Traits::FluxVariablesCache;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = false;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    FaceCenteredStaggeredGridFluxVariablesCache(const Problem& problem) : problemPtr_(&problem) {}

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
