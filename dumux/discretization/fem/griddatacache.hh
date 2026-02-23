// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FEDiscretization
 * \brief Global data cache
 */
#ifndef DUMUX_DISCRETIZATION_FE_GRID_DATACACHE_HH
#define DUMUX_DISCRETIZATION_FE_GRID_DATACACHE_HH

#include <memory>
#include <ranges>
#include <unordered_map>

#include <dumux/parallel/parallel_for.hh>

// make the local view function available whenever we use this class
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/fem/elementdatacache.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

namespace Dumux {

/*!
 * \ingroup FEDiscretization
 * \brief Data cache traits
 */
template<class P, class DC>
struct FEDefaultGridDataCacheTraits
{
    using Problem = P;
    using DataCache = DC;

    template<class GridDataCache, bool cachingEnabled>
    using LocalView = FEElementDataCache<GridDataCache, cachingEnabled>;
};

/*!
 * \ingroup FEDiscretization
 * \brief Data cache implementation on a gridview
 * \note The class is specialized for a version with and without grid caching, and for different quadrature rules
 */
template<class Problem,
         class DataCache,
         bool cachingEnabled = false,
         class Traits = FEDefaultGridDataCacheTraits<Problem, DataCache>>
class FEGridDataCache;

/*!
 * \ingroup FEDiscretization
 * \brief Data caches on a gridview with grid caching enabled (general quadrature specialization)
 * \note The data caches of the gridview are stored which is memory intensive but faster
 */
template<class P, class DC, class Traits>
class FEGridDataCache<P, DC, true, Traits>
{
    using Problem = typename Traits::Problem;
    using ThisType = FEGridDataCache<P, DC, true, Traits>;

public:
    //! export the local data cache type
    using DataCache = typename Traits::DataCache;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = true;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    FEGridDataCache(const Problem& problem) : problemPtr_(&problem) {}

    template<class GridGeometry, class GridVariables, class SolutionVector>
    void update(const GridGeometry& gridGeometry,
                const GridVariables& gridVars,
                const SolutionVector& sol,
                bool forceUpdate = false)
    {
        // Here, we do not do anything unless it is a forced update
        if (forceUpdate)
        {
            elementCache_.resize(gridGeometry.gridView().size(0));
            boundaryIntersectionCache_.resize(gridGeometry.gridView().size(0));
            Dumux::parallelFor(gridGeometry.gridView().size(0), [&, &problem = problem()](const std::size_t eIdx)
            {
                // Prepare the geometries within the elements of the stencil
                const auto element = gridGeometry.element(eIdx);
                const auto eGeometryView = localView(gridGeometry).bind(element);
                const auto elemVars = localView(gridVars).bind(element, eGeometryView, sol);

                // Resize element cache based on element quadrature rule
                const auto elemQuadRule = CVFE::quadratureRule(eGeometryView, element);
                elementCache_[eIdx].resize(std::ranges::size(elemQuadRule));
                for (const auto& qpData : elemQuadRule)
                    elementCache_[eIdx][qpData.ipData().qpIndex()].update(problem, element, eGeometryView, elemVars, qpData.ipData());

                // Rebuild boundary intersection cache for this element
                if (!boundaryIntersectionCache_[eIdx])
                    boundaryIntersectionCache_[eIdx] = std::make_unique<std::unordered_map<int, std::vector<DataCache>>>();
                else
                    boundaryIntersectionCache_[eIdx]->clear();

                for (const auto& intersection : intersections(gridGeometry.gridView(), element))
                {
                    if (intersection.boundary())
                    {
                        const auto quadRule = CVFE::quadratureRule(eGeometryView, intersection);
                        const auto iIdx = intersection.indexInInside();
                        (*boundaryIntersectionCache_[eIdx])[iIdx].resize(std::ranges::size(quadRule));

                        for (const auto& qpData : quadRule)
                            (*boundaryIntersectionCache_[eIdx])[iIdx][qpData.ipData().qpIndex()].update(problem,
                                                                                                        element,
                                                                                                        eGeometryView,
                                                                                                        elemVars,
                                                                                                        qpData.ipData());
                    }
                }
            });
        }
    }

    template<class ElementGeometryView, class ElementVariables>
    void updateElement(const typename ElementGeometryView::Element& element,
                       const ElementGeometryView& eGeometryView,
                       const ElementVariables& elemVars)
    {
        if constexpr (DataCache::isSolDependent)
        {
            const auto eIdx = eGeometryView.gridGeometry().elementMapper().index(element);

            // Resize element cache based on element quadrature rule
            const auto elemQuadRule = CVFE::quadratureRule(eGeometryView, element);
            elementCache_[eIdx].resize(std::ranges::size(elemQuadRule));
            for (const auto& qpData : elemQuadRule)
                elementCache_[eIdx][qpData.ipData().qpIndex()].update(problem(), element, eGeometryView, elemVars, qpData.ipData());

            // Rebuild boundary intersection cache for this element
            if (!boundaryIntersectionCache_[eIdx])
                boundaryIntersectionCache_[eIdx] = std::make_unique<std::unordered_map<int, std::vector<DataCache>>>();
            else
                boundaryIntersectionCache_[eIdx]->clear();

            for (const auto& intersection : intersections(eGeometryView.gridGeometry().gridView(), element))
            {
                if (intersection.boundary())
                {
                    const auto quadRule = CVFE::quadratureRule(eGeometryView, intersection);
                    const auto iIdx = intersection.indexInInside();
                    (*boundaryIntersectionCache_[eIdx])[iIdx].resize(std::ranges::size(quadRule));

                    for (const auto& qpData : quadRule)
                        (*boundaryIntersectionCache_[eIdx])[iIdx][qpData.ipData().qpIndex()].update(problem(),
                                                                                                    element,
                                                                                                    eGeometryView,
                                                                                                    elemVars,
                                                                                                    qpData.ipData());
                }
            }
        }
    }

    const Problem& problem() const
    { return *problemPtr_; }

    // access operator
    const DataCache& elementCache(std::size_t eIdx, std::size_t qpIdx) const { return elementCache_[eIdx][qpIdx]; }
    DataCache& elementCache(std::size_t eIdx, std::size_t qpIdx) { return elementCache_[eIdx][qpIdx]; }

    // access operator for boundary intersection cache
    DataCache& boundaryIntersectionCache(std::size_t eIdx, int iIdx, std::size_t qpIdx) { return (*boundaryIntersectionCache_[eIdx])[iIdx][qpIdx]; }
    const DataCache& boundaryIntersectionCache(std::size_t eIdx, int iIdx, std::size_t qpIdx) const { return (*boundaryIntersectionCache_[eIdx]).at(iIdx)[qpIdx]; }

private:
    // currently bound element
    const Problem* problemPtr_;
    std::vector<std::vector<DataCache>> elementCache_; //! storage per element/qp
    std::vector<std::unique_ptr<std::unordered_map<int, std::vector<DataCache>>>> boundaryIntersectionCache_; //! storage per element/boundary intersection/qp
};

/*!
 * \ingroup FEDiscretization
 * \brief Data caches on a gridview with grid caching disabled
 */
template<class P, class DC, class Traits>
class FEGridDataCache<P, DC, false, Traits>
{
    using Problem = typename Traits::Problem;
    using ThisType = FEGridDataCache<P, DC, false, Traits>;

public:
    //! export the local data cache type
    using DataCache = typename Traits::DataCache;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = false;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    FEGridDataCache(const Problem& problem) : problemPtr_(&problem) {}

    template<class GridGeometry, class GridVariables, class SolutionVector>
    void update(const GridGeometry& gridGeometry,
                const GridVariables& gridVars,
                const SolutionVector& sol,
                bool forceUpdate = false) {}

    const Problem& problem() const
    { return *problemPtr_; }

private:
    const Problem* problemPtr_;
};

} // end namespace Dumux

#endif
