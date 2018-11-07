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
 * \ingroup GodunovDiscretization
 * \brief Flux variable caches on a gridview
 */
#ifndef DUMUX_DISCRETIZATION_GODUNOV_GRID_FLUXVARSCACHE_HH
#define DUMUX_DISCRETIZATION_GODUNOV_GRID_FLUXVARSCACHE_HH

//! make the local view function available whenever we use this class
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/cellcentered/godunov/elementfluxvariablescache.hh>

namespace Dumux {

/*!
 * \ingroup GodunovDiscretization
 * \brief Flux variable caches traits
 */
template<class P, class FVC, class FVCF>
struct GodunovDefaultGridFVCTraits
{
    using Problem = P;
    using FluxVariablesCache = FVC;
    using FluxVariablesCacheFiller = FVCF;

    template<class GridFluxVariablesCache, bool cachingEnabled>
    using LocalView = GodunovElementFluxVariablesCache<GridFluxVariablesCache, cachingEnabled>;
};

/*!
 * \ingroup GodunovDiscretization
 * \brief Flux variable caches on a gridview
 * \note The class is specialized for a version with and without grid caching
 */
template<class Problem,
         class FluxVariablesCache,
         class FluxVariablesCacheFiller,
         bool EnableGridFluxVariablesCache = false,
         class Traits = GodunovDefaultGridFVCTraits<Problem, FluxVariablesCache, FluxVariablesCacheFiller> >
class GodunovGridFluxVariablesCache;

/*!
 * \ingroup GodunovDiscretization
 * \brief Flux variable caches on a gridview with grid caching enabled
 * \note The flux caches of the gridview are stored which is memory intensive but faster
 */
template<class P, class FVC, class FVCF, class Traits>
class GodunovGridFluxVariablesCache<P, FVC, FVCF, true, Traits>
{
    using Problem = typename Traits::Problem;
    using ThisType = GodunovGridFluxVariablesCache<P, FVC, FVCF, true, Traits>;

public:
    //! export the flux variable cache type
    using FluxVariablesCache = typename Traits::FluxVariablesCache;

    //! export the flux variable cache filler type
    using FluxVariablesCacheFiller = typename Traits::FluxVariablesCacheFiller;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = true;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    // The constructor
    GodunovGridFluxVariablesCache(const Problem& problem) : problemPtr_(&problem) {}

    // When global caching is enabled, precompute transmissibilities and stencils for all the scv faces
    template<class FVGridGeometry, class GridVolumeVariables, class SolutionVector>
    void update(const FVGridGeometry& fvGridGeometry,
                const GridVolumeVariables& gridVolVars,
                const SolutionVector& sol,
                bool forceUpdate = false)
    {
        // only do the update if fluxes are solution dependent or if update is forced
        if (FluxVariablesCacheFiller::isSolDependent || forceUpdate)
        {
            // instantiate helper class to fill the caches
            FluxVariablesCacheFiller filler(problem());

            fluxVarsCache_.resize(fvGridGeometry.numScvf());
            for (const auto& element : elements(fvGridGeometry.gridView()))
            {
                // Prepare the geometries within the elements of the stencil
                auto fvGeometry = localView(fvGridGeometry);
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

    template<class FVElementGeometry, class ElementVolumeVariables>
    void updateElement(const typename FVElementGeometry::FVGridGeometry::GridView::template Codim<0>::Entity& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars)
    {
        if (FluxVariablesCacheFiller::isSolDependent)
        {
            const auto globalI = fvGeometry.fvGridGeometry().elementMapper().index(element);

            // instantiate filler class
            FluxVariablesCacheFiller filler(problem());

            // update the caches inside this element
            for (const auto& scvf : scvfs(fvGeometry))
                filler.fill(*this, fluxVarsCache_[scvf.index()], element, fvGeometry, elemVolVars, scvf);

            // update the caches in the neighbors
            for (const auto& dataJ : fvGeometry.fvGridGeometry().connectivityMap()[globalI])
            {
                const auto elementJ = fvGeometry.fvGridGeometry().element(dataJ.globalJ);
                for (const auto scvfIdxJ : dataJ.scvfsJ)
                {
                    const auto& scvfJ = fvGeometry.scvf(scvfIdxJ);
                    filler.fill(*this, fluxVarsCache_[scvfJ.index()], elementJ, fvGeometry, elemVolVars, scvfJ);
                }
            }
        }
    }

    // access operators in the case of caching
    template<class SubControlVolumeFace>
    const FluxVariablesCache& operator [](const SubControlVolumeFace& scvf) const
    { return fluxVarsCache_[scvf.index()]; }

    template<class SubControlVolumeFace>
    FluxVariablesCache& operator [](const SubControlVolumeFace& scvf)
    { return fluxVarsCache_[scvf.index()]; }

    const Problem& problem() const
    { return *problemPtr_; }

private:
    const Problem* problemPtr_;

    std::vector<FluxVariablesCache> fluxVarsCache_;
    std::vector<std::size_t> globalScvfIndices_;
};

/*!
 * \ingroup GodunovDiscretization
 * \brief Flux variable caches on a gridview with grid caching disabled
 */
template<class P, class FVC, class FVCF, class Traits>
class GodunovGridFluxVariablesCache<P, FVC, FVCF, false, Traits>
{
    using Problem = typename Traits::Problem;
    using ThisType = GodunovGridFluxVariablesCache<P, FVC, FVCF, false, Traits>;

public:
    //! export the flux variable cache type
    using FluxVariablesCache = typename Traits::FluxVariablesCache;

    //! export the flux variable cache filler type
    using FluxVariablesCacheFiller = typename Traits::FluxVariablesCacheFiller;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = false;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    // The constructor
    GodunovGridFluxVariablesCache(const Problem& problem) : problemPtr_(&problem) {}

    //! When global flux variables caching is disabled, we don't need to update the cache
    template<class FVGridGeometry, class GridVolumeVariables, class SolutionVector>
    void update(const FVGridGeometry& fvGridGeometry,
                const GridVolumeVariables& gridVolVars,
                const SolutionVector& sol,
                bool forceUpdate = false) {}

    //! When global flux variables caching is disabled, we don't need to update the cache
    template<class FVElementGeometry, class ElementVolumeVariables>
    void updateElement(const typename FVElementGeometry::FVGridGeometry::GridView::template Codim<0>::Entity& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars) {}

    const Problem& problem() const
    { return *problemPtr_; }

private:
    const Problem* problemPtr_;
};

} // end namespace Dumux

#endif
