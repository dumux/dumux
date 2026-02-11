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
#ifndef DUMUX_DISCRETIZATION_HYBRID_CVFE_GRID_FLUXVARSCACHE_HH
#define DUMUX_DISCRETIZATION_HYBRID_CVFE_GRID_FLUXVARSCACHE_HH

#include <memory>
#include <ranges>
#include <unordered_map>

#include <dumux/parallel/parallel_for.hh>

// make the local view function available whenever we use this class
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/cvfe/hybrid/elementfluxvariablescache.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

namespace Dumux {

/*!
 * \ingroup CVFEDiscretization
 * \brief Flux variable caches traits
 */
template<class P, class FVC>
struct HybridCVFEDefaultGridFVCTraits
{
    using Problem = P;
    using FluxVariablesCache = FVC;

    template<class GridFluxVariablesCache, bool cachingEnabled>
    using LocalView = HybridCVFEElementFluxVariablesCache<GridFluxVariablesCache, cachingEnabled>;
};

/*!
 * \ingroup CVFEDiscretization
 * \brief Flux variable caches implementation on a gridview
 * \note The class is specialized for a version with and without grid caching, and for different quadrature rules
 */
template<class Problem,
         class FluxVariablesCache,
         bool cachingEnabled = false,
         class Traits = HybridCVFEDefaultGridFVCTraits<Problem, FluxVariablesCache>>
class HybridCVFEGridFluxVariablesCache;

/*!
 * \ingroup CVFEDiscretization
 * \brief Flux variable caches on a gridview with grid caching enabled (general quadrature specialization)
 * \note The flux caches of the gridview are stored which is memory intensive but faster
 * \note This implementation supports multiple quadrature points per sub-control volume face
 */
template<class P, class FVC, class Traits>
class HybridCVFEGridFluxVariablesCache<P, FVC, true, Traits>
{
    using Problem = typename Traits::Problem;
    using ThisType = HybridCVFEGridFluxVariablesCache<P, FVC, true, Traits>;

public:
    //! export the flux variable cache type
    using FluxVariablesCache = typename Traits::FluxVariablesCache;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = true;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    HybridCVFEGridFluxVariablesCache(const Problem& problem) : problemPtr_(&problem) {}

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
            scvfOffset_.resize(gridGeometry.gridView().size(0));
            elementCache_.resize(gridGeometry.gridView().size(0));
            boundaryIntersectionCache_.resize(gridGeometry.gridView().size(0));
            Dumux::parallelFor(gridGeometry.gridView().size(0), [&, &problem = problem()](const std::size_t eIdx)
            {
                // Prepare the geometries within the elements of the stencil
                const auto element = gridGeometry.element(eIdx);
                const auto fvGeometry = localView(gridGeometry).bind(element);
                const auto elemVolVars = localView(gridVolVars).bind(element, fvGeometry, sol);

                // Build offset vector and resize flat cache
                scvfOffset_[eIdx].resize(fvGeometry.numScvf() + 1, 0);
                for (const auto& scvf : scvfs(fvGeometry))
                {
                    const auto numQps = std::ranges::size(CVFE::quadratureRule(fvGeometry, scvf));
                    scvfOffset_[eIdx][scvf.index() + 1] = numQps;
                }
                for (std::size_t i = 2; i < scvfOffset_[eIdx].size(); ++i)
                    scvfOffset_[eIdx][i] += scvfOffset_[eIdx][i-1];

                fluxVarsCache_[eIdx].resize(scvfOffset_[eIdx].back());
                for (const auto& scvf : scvfs(fvGeometry))
                {
                    for (const auto& qpData : CVFE::quadratureRule(fvGeometry, scvf))
                        cache(eIdx, qpData.ipData().scvfIndex(), qpData.ipData().qpIndex()).update(problem,
                                                                                                   element,
                                                                                                   fvGeometry,
                                                                                                   elemVolVars,
                                                                                                   qpData.ipData());
                }

                // Resize element cache based on element quadrature rule
                const auto elemQuadRule = CVFE::quadratureRule(fvGeometry, element);
                elementCache_[eIdx].resize(std::ranges::size(elemQuadRule));
                for (const auto& qpData : elemQuadRule)
                    elementCache_[eIdx][qpData.ipData().qpIndex()].update(problem, element, fvGeometry, elemVolVars, qpData.ipData());

                // Rebuild boundary intersection cache for this element
                if (!boundaryIntersectionCache_[eIdx])
                    boundaryIntersectionCache_[eIdx] = std::make_unique<std::unordered_map<int, std::vector<FluxVariablesCache>>>();
                else
                    boundaryIntersectionCache_[eIdx]->clear();

                for (const auto& intersection : intersections(gridGeometry.gridView(), element))
                {
                    if (intersection.boundary())
                    {
                        const auto quadRule = CVFE::quadratureRule(fvGeometry, intersection);
                        const auto iIdx = intersection.indexInInside();
                        (*boundaryIntersectionCache_[eIdx])[iIdx].resize(std::ranges::size(quadRule));

                        for (const auto& qpData : quadRule)
                            (*boundaryIntersectionCache_[eIdx])[iIdx][qpData.ipData().qpIndex()].update(problem,
                                                                                                        element,
                                                                                                        fvGeometry,
                                                                                                        elemVolVars,
                                                                                                        qpData.ipData());
                    }
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

            // Build offset vector and resize flat cache
            scvfOffset_[eIdx].resize(fvGeometry.numScvf() + 1, 0);
            for (const auto& scvf : scvfs(fvGeometry))
            {
                const auto numQps = std::ranges::size(CVFE::quadratureRule(fvGeometry, scvf));
                scvfOffset_[eIdx][scvf.index() + 1] = numQps;
            }
            for (std::size_t i = 2; i < scvfOffset_[eIdx].size(); ++i)
                scvfOffset_[eIdx][i] += scvfOffset_[eIdx][i-1];

            fluxVarsCache_[eIdx].resize(scvfOffset_[eIdx].back());
            for (const auto& scvf : scvfs(fvGeometry))
            {
                for (const auto& qpData : CVFE::quadratureRule(fvGeometry, scvf))
                    cache(eIdx, qpData.ipData().scvfIndex(), qpData.ipData().qpIndex()).update(problem(),
                                                                                               element,
                                                                                               fvGeometry,
                                                                                               elemVolVars,
                                                                                               qpData.ipData());
            }

            // Resize element cache based on element quadrature rule
            const auto elemQuadRule = CVFE::quadratureRule(fvGeometry, element);
            elementCache_[eIdx].resize(std::ranges::size(elemQuadRule));
            for (const auto& qpData : elemQuadRule)
                elementCache_[eIdx][qpData.ipData().qpIndex()].update(problem(), element, fvGeometry, elemVolVars, qpData.ipData());

            // Rebuild boundary intersection cache for this element
            if (!boundaryIntersectionCache_[eIdx])
                boundaryIntersectionCache_[eIdx] = std::make_unique<std::unordered_map<int, std::vector<FluxVariablesCache>>>();
            else
                boundaryIntersectionCache_[eIdx]->clear();

            for (const auto& intersection : intersections(fvGeometry.gridGeometry().gridView(), element))
            {
                if (intersection.boundary())
                {
                    const auto quadRule = CVFE::quadratureRule(fvGeometry, intersection);
                    const auto iIdx = intersection.indexInInside();
                    (*boundaryIntersectionCache_[eIdx])[iIdx].resize(std::ranges::size(quadRule));

                    for (const auto& qpData : quadRule)
                        (*boundaryIntersectionCache_[eIdx])[iIdx][qpData.ipData().qpIndex()].update(problem(),
                                                                                                    element,
                                                                                                    fvGeometry,
                                                                                                    elemVolVars,
                                                                                                    qpData.ipData());
                }
            }
        }
    }

    const Problem& problem() const
    { return *problemPtr_; }

    // access operator
    const FluxVariablesCache& cache(std::size_t eIdx, std::size_t scvfIdx, std::size_t qpIdx) const
    { return fluxVarsCache_[eIdx][scvfOffset_[eIdx][scvfIdx] + qpIdx]; }

    // access operator
    FluxVariablesCache& cache(std::size_t eIdx, std::size_t scvfIdx, std::size_t qpIdx)
    { return fluxVarsCache_[eIdx][scvfOffset_[eIdx][scvfIdx] + qpIdx]; }

    // access operator
    const FluxVariablesCache& elementCache(std::size_t eIdx, std::size_t qpIdx) const { return elementCache_[eIdx][qpIdx]; }
    FluxVariablesCache& elementCache(std::size_t eIdx, std::size_t qpIdx) { return elementCache_[eIdx][qpIdx]; }

    // access operator for boundary intersection cache
    FluxVariablesCache& boundaryIntersectionCache(std::size_t eIdx, int iIdx, std::size_t qpIdx) { return (*boundaryIntersectionCache_[eIdx])[iIdx][qpIdx]; }
    const FluxVariablesCache& boundaryIntersectionCache(std::size_t eIdx, int iIdx, std::size_t qpIdx) const { return (*boundaryIntersectionCache_[eIdx]).at(iIdx)[qpIdx]; }

private:
    // currently bound element
    const Problem* problemPtr_;
    std::vector<std::vector<FluxVariablesCache>> fluxVarsCache_; //! flat storage per element
    std::vector<std::vector<std::size_t>> scvfOffset_; //! scvfOffset_[eIdx][scvfIdx] -> offset in fluxVarsCache_[eIdx]
    std::vector<std::vector<FluxVariablesCache>> elementCache_; //! storage per element/qp
    std::vector<std::unique_ptr<std::unordered_map<int, std::vector<FluxVariablesCache>>>> boundaryIntersectionCache_; //! storage per element/boundary intersection/qp
};

/*!
 * \ingroup CVFEDiscretization
 * \brief Flux variable caches on a gridview with grid caching disabled
 */
template<class P, class FVC, class Traits>
class HybridCVFEGridFluxVariablesCache<P, FVC, false, Traits>
{
    using Problem = typename Traits::Problem;
    using ThisType = HybridCVFEGridFluxVariablesCache<P, FVC, false, Traits>;

public:
    //! export the flux variable cache type
    using FluxVariablesCache = typename Traits::FluxVariablesCache;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = false;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    HybridCVFEGridFluxVariablesCache(const Problem& problem) : problemPtr_(&problem) {}

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
