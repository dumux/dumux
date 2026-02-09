// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
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
#include <dumux/discretization/cvfe/quadraturerules.hh>

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
 * \brief Flux variable caches implementation on a gridview
 * \note The class is specialized for a version with and without grid caching, and for different quadrature rules
 */
template<class Problem,
         class FluxVariablesCache,
         bool cachingEnabled,
         class Traits,
         class ScvfQuadratureRule>
class CVFEGridFluxVariablesCacheImpl;

/*!
 * \ingroup CVFEDiscretization
 * \brief Flux variable caches on a gridview
 */
template<class Problem,
         class FluxVariablesCache,
         bool cachingEnabled = false,
         class Traits = CVFEDefaultGridFVCTraits<Problem, FluxVariablesCache>>
using CVFEGridFluxVariablesCache = CVFEGridFluxVariablesCacheImpl<Problem, FluxVariablesCache, cachingEnabled,
                                                                  Traits, Detail::ScvfQuadratureRuleOrDefault_t<typename Traits::FluxVariablesCache>>;

/*!
 * \ingroup CVFEDiscretization
 * \brief Flux variable caches on a gridview with grid caching enabled (MidpointQuadrature specialization)
 * \note The flux caches of the gridview are stored which is memory intensive but faster
 * \note This specialization uses the implementation for midpoint quadrature rule
 */
template<class P, class FVC, class Traits>
class CVFEGridFluxVariablesCacheImpl<P, FVC, true, Traits, QuadratureRules::MidpointQuadrature>
{
    using Problem = typename Traits::Problem;
    using ThisType = CVFEGridFluxVariablesCacheImpl<P, FVC, true, Traits, QuadratureRules::MidpointQuadrature>;

public:
    //! export the flux variable cache type
    using FluxVariablesCache = typename Traits::FluxVariablesCache;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = true;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    CVFEGridFluxVariablesCacheImpl(const Problem& problem) : problemPtr_(&problem) {}

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
 * \brief Flux variable caches on a gridview with grid caching enabled (general quadrature specialization)
 * \note The flux caches of the gridview are stored which is memory intensive but faster
 * \note This implementation supports multiple quadrature points per sub-control volume face
 */
template<class P, class FVC, class Traits, class ScvfQR>
class CVFEGridFluxVariablesCacheImpl<P, FVC, true, Traits, ScvfQR>
{
    using Problem = typename Traits::Problem;
    using ThisType = CVFEGridFluxVariablesCacheImpl<P, FVC, true, Traits, ScvfQR>;

public:
    //! export the flux variable cache type
    using FluxVariablesCache = typename Traits::FluxVariablesCache;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = true;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    CVFEGridFluxVariablesCacheImpl(const Problem& problem) : problemPtr_(&problem) {}

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
            qpsOffset_.resize(gridGeometry.gridView().size(0));
            Dumux::parallelFor(gridGeometry.gridView().size(0), [&, &problem = problem()](const std::size_t eIdx)
            {
                // Prepare the geometries within the elements of the stencil
                const auto element = gridGeometry.element(eIdx);
                const auto fvGeometry = localView(gridGeometry).bind(element);
                const auto elemVolVars = localView(gridVolVars).bind(element, fvGeometry, sol);

                // Compute total number of cache entries and offsets for this element
                qpsOffset_[eIdx].resize(fvGeometry.numScvf() + 1, 0);
                for (const auto& scvf : scvfs(fvGeometry))
                {
                    const auto numQps = std::ranges::size(CVFE::quadratureRule(fvGeometry, scvf));
                    qpsOffset_[eIdx][scvf.index() + 1] = numQps;
                }
                for (std::size_t i = 2; i < qpsOffset_[eIdx].size(); ++i)
                    qpsOffset_[eIdx][i] += qpsOffset_[eIdx][i-1];

                // Resize flat cache and update entries
                fluxVarsCache_[eIdx].resize(qpsOffset_[eIdx].back());
                for (const auto& scvf : scvfs(fvGeometry))
                {
                    for (const auto& qpData : CVFE::quadratureRule(fvGeometry, scvf))
                        cache(eIdx, qpData.ipData().scvfIndex(), qpData.ipData().qpIndex()).update(problem,
                                                                                                   element,
                                                                                                   fvGeometry,
                                                                                                   elemVolVars,
                                                                                                   qpData.ipData());
                }
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

            // Compute offsets for this element
            qpsOffset_[eIdx].resize(fvGeometry.numScvf() + 1, 0);
            for (const auto& scvf : scvfs(fvGeometry))
            {
                const auto numQps = std::ranges::size(CVFE::quadratureRule(fvGeometry, scvf));
                qpsOffset_[eIdx][scvf.index() + 1] = numQps;
            }
            for (std::size_t i = 2; i < qpsOffset_[eIdx].size(); ++i)
                qpsOffset_[eIdx][i] += qpsOffset_[eIdx][i-1];

            // Resize and update flat cache
            fluxVarsCache_[eIdx].resize(qpsOffset_[eIdx].back());
            for (const auto& scvf : scvfs(fvGeometry))
            {
                for (const auto& qpData : CVFE::quadratureRule(fvGeometry, scvf))
                    cache(eIdx, qpData.ipData().scvfIndex(), qpData.ipData().qpIndex()).update(problem(),
                                                                                               element,
                                                                                               fvGeometry,
                                                                                               elemVolVars,
                                                                                               qpData.ipData());
            }
        }
    }

    const Problem& problem() const
    { return *problemPtr_; }

    // access operator
    const FluxVariablesCache& cache(std::size_t eIdx, std::size_t scvfIdx, std::size_t qpIdx) const
    { return fluxVarsCache_[eIdx][qpsOffset_[eIdx][scvfIdx] + qpIdx]; }

    // access operator
    FluxVariablesCache& cache(std::size_t eIdx, std::size_t scvfIdx, std::size_t qpIdx)
    { return fluxVarsCache_[eIdx][qpsOffset_[eIdx][scvfIdx] + qpIdx]; }

private:
    // currently bound element
    const Problem* problemPtr_;
    std::vector<std::vector<FluxVariablesCache>> fluxVarsCache_; //! storage per element
    std::vector<std::vector<std::size_t>> qpsOffset_; //! offset[eIdx][scvfIdx] -> start index in fluxVarsCache_[eIdx]
};

/*!
 * \ingroup CVFEDiscretization
 * \brief Flux variable caches on a gridview with grid caching disabled
 */
template<class P, class FVC, class Traits, class ScvfQR>
class CVFEGridFluxVariablesCacheImpl<P, FVC, false, Traits, ScvfQR>
{
    using Problem = typename Traits::Problem;
    using ThisType = CVFEGridFluxVariablesCacheImpl<P, FVC, false, Traits, ScvfQR>;

public:
    //! export the flux variable cache type
    using FluxVariablesCache = typename Traits::FluxVariablesCache;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = false;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    CVFEGridFluxVariablesCacheImpl(const Problem& problem) : problemPtr_(&problem) {}

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
