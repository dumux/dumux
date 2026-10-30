// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CVFEDiscretization
 * \brief The grid variables cache class for control-volume finite element methods
 */
#ifndef DUMUX_DISCRETIZATION_CVFE_GRID_VARIABLES_CACHE_HH
#define DUMUX_DISCRETIZATION_CVFE_GRID_VARIABLES_CACHE_HH

#include <vector>
#include <ranges>
#include <memory>

#include <dumux/parallel/parallel_for.hh>

#include <dumux/common/concepts/localdofs_.hh>

// make the local view function available whenever we use this class
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/cvfe/elementvariables_.hh>
#include <dumux/discretization/cvfe/elementsolution.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

namespace Dumux::Experimental::CVFE {

template<class P, class V, class IPD>
struct CVFEDefaultGridVariablesCacheTraits
{
    using Problem = P;
    using Variables = V;
    using InterpolationPointData = IPD;

    template<class GridVariablesCache, bool cachingEnabled>
    using LocalView = CVFEElementVariables<GridVariablesCache, cachingEnabled>;
};

/*!
 * \ingroup CVFEDiscretization
 * \brief The grid variables cache class for control-volume finite element methods
 * \tparam Traits the traits
 * \tparam enableCaching whether the variables are stored for the whole grid
 * \tparam enableIpDataCaching whether the interpolation point data are stored for the whole grid
 */
template<class Traits, bool enableCaching, bool enableIpDataCaching = enableCaching>
class CVFEGridVariablesCache;

// specialization in case of storing the local variables
template<class Traits, bool enableIpDataCaching>
class CVFEGridVariablesCache<Traits, /*cachingEnabled*/true, enableIpDataCaching>
{
    using ThisType = CVFEGridVariablesCache<Traits, true, enableIpDataCaching>;
public:
    //! export the problem type
    using Problem = typename Traits::Problem;

    //! export the variables type
    using Variables = typename Traits::Variables;

    //! export interpolation point data type
    using InterpolationPointData = typename Traits::InterpolationPointData;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = true;

    //! make it possible to query if the interpolation point data are stored for the whole grid
    static constexpr bool interpolationPointDataCachingEnabled = enableIpDataCaching;

    static_assert(enableIpDataCaching || !InterpolationPointData::isSolDependent,
                  "Solution-dependent interpolation point data are stored together with the variables");

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    //! export the type of the mutable local view
    using MutableLocalView = typename LocalView::MutableView;

    CVFEGridVariablesCache(const Problem& problem) : problemPtr_(&problem) {}

    template<class GridDiscretization, class SolutionVector>
    void init(const GridDiscretization& gridDiscretization, const SolutionVector& sol)
    {
        variables_.resize(gridDiscretization.gridView().size(0));
        if constexpr (enableIpDataCaching)
            ipDataCache_ = std::make_shared<std::vector<ElementInterpolationPointDataCache>>(gridDiscretization.gridView().size(0));

        Dumux::parallelFor(gridDiscretization.gridView().size(0), [&, &problem = problem()](const std::size_t eIdx)
        {
            const auto element = gridDiscretization.element(eIdx);
            const auto elemDisc = localView(gridDiscretization).bindElement(element);

            // get the element solution
            auto elemSol = elementSolution(element, sol, gridDiscretization);

            variables_[eIdx].resize(Dumux::Detail::LocalDofs::numLocalDofs(elemDisc));
            for (const auto& localDof : localDofs(elemDisc))
                variables_[eIdx][localDof.index()].update(elemSol, problem, elemDisc, ipData(elemDisc, localDof));

            if constexpr (enableIpDataCaching)
                (*ipDataCache_)[eIdx].update(problem, element, elemDisc, variables_[eIdx]);
        });
    }

    template<class GridDiscretization, class SolutionVector>
    void update(const GridDiscretization& gridDiscretization, const SolutionVector& sol)
    {
        if constexpr (InterpolationPointData::isSolDependent)
        {
            auto newIpDataCache = std::make_shared<std::vector<ElementInterpolationPointDataCache>>(*ipDataCache_);

            Dumux::parallelFor(gridDiscretization.gridView().size(0), [&, &problem = problem(), newIpDataCache](const std::size_t eIdx)
            {
                const auto element = gridDiscretization.element(eIdx);
                const auto elemDisc = localView(gridDiscretization).bindElement(element);

                // get the element solution
                auto elemSol = elementSolution(element, sol, gridDiscretization);

                for (const auto& localDof : localDofs(elemDisc))
                    variables_[eIdx][localDof.index()].update(elemSol, problem, elemDisc, ipData(elemDisc, localDof));

                (*newIpDataCache)[eIdx].update(problem, element, elemDisc, variables_[eIdx]);
            });

            ipDataCache_ = std::move(newIpDataCache);
        }
        else
        {
            Dumux::parallelFor(gridDiscretization.gridView().size(0), [&, &problem = problem()](const std::size_t eIdx)
            {
                const auto element = gridDiscretization.element(eIdx);
                const auto elemDisc = localView(gridDiscretization).bindElement(element);

                // get the element solution
                auto elemSol = elementSolution(element, sol, gridDiscretization);

                for (const auto& localDof : localDofs(elemDisc))
                    variables_[eIdx][localDof.index()].update(elemSol, problem, elemDisc, ipData(elemDisc, localDof));
            });
        }
    }

    template<class ScvOrLocalDof>
    const Variables& variables(const ScvOrLocalDof& scvOrLocalDof) const
    {
        if constexpr (Concept::LocalDof<ScvOrLocalDof>)
            return variables_[scvOrLocalDof.elementIndex()][scvOrLocalDof.index()];
        else
            return variables_[scvOrLocalDof.elementIndex()][scvOrLocalDof.localDofIndex()];
    }

    template<class ScvOrLocalDof>
    Variables& variables(const ScvOrLocalDof& scvOrLocalDof)
    {
        if constexpr (Concept::LocalDof<ScvOrLocalDof>)
            return variables_[scvOrLocalDof.elementIndex()][scvOrLocalDof.index()];
        else
            return variables_[scvOrLocalDof.elementIndex()][scvOrLocalDof.localDofIndex()];
    }

    const Variables& variables(const std::size_t eIdx, const std::size_t localIdx) const
    { return variables_[eIdx][localIdx]; }

    Variables& variables(const std::size_t eIdx, const std::size_t localIdx)
    { return variables_[eIdx][localIdx]; }

    //! The variables of all local dofs of an element
    const std::vector<Variables>& localDofVariables(const std::size_t eIdx) const
    { return variables_[eIdx]; }

    const InterpolationPointData& scvfCache(std::size_t eIdx, std::size_t scvfIdx, std::size_t qpIdx) const
    { return (*ipDataCache_)[eIdx].scvfCache(scvfIdx, qpIdx); }

    InterpolationPointData& scvfCache(std::size_t eIdx, std::size_t scvfIdx, std::size_t qpIdx)
    { return (*ipDataCache_)[eIdx].scvfCache(scvfIdx, qpIdx); }

    const auto& cache(std::size_t eIdx) const
    { return (*ipDataCache_)[eIdx]; }

    auto& cache(std::size_t eIdx)
    { return (*ipDataCache_)[eIdx]; }

    const Problem& problem() const
    { return *problemPtr_; }

private:
    using ElementInterpolationPointDataCache = Dumux::Detail::CVFE::ScvfInterpolationPointDataCache<InterpolationPointData>;

    const Problem* problemPtr_;
    std::vector<std::vector<Variables>> variables_;
    std::shared_ptr<std::vector<ElementInterpolationPointDataCache>> ipDataCache_;
};

// Specialization when the current local variables are not stored
template<class Traits, bool enableIpDataCaching>
class CVFEGridVariablesCache<Traits, /*cachingEnabled*/false, enableIpDataCaching>
{
    using ThisType = CVFEGridVariablesCache<Traits, false, enableIpDataCaching>;

public:
    //! export the problem type
    using Problem = typename Traits::Problem;

    //! export the variables type
    using Variables = typename Traits::Variables;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = false;

    //! export interpolation point data type
    using InterpolationPointData = typename Traits::InterpolationPointData;

    //! make it possible to query if the interpolation point data are stored for the whole grid
    static constexpr bool interpolationPointDataCachingEnabled = enableIpDataCaching;

    static_assert(!enableIpDataCaching || !InterpolationPointData::isSolDependent,
                  "Solution-dependent interpolation point data are stored together with the variables");

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    //! export the type of the mutable local view
    using MutableLocalView = typename LocalView::MutableView;

    CVFEGridVariablesCache(const Problem& problem) : problemPtr_(&problem) {}

    template<class GridDiscretization, class SolutionVector>
    void init(const GridDiscretization& gridDiscretization, const SolutionVector& sol)
    {
        if constexpr (enableIpDataCaching)
        {
            auto ipDataCache = std::make_shared<std::vector<ElementInterpolationPointDataCache>>(gridDiscretization.gridView().size(0));
            Dumux::parallelFor(gridDiscretization.gridView().size(0), [&, &problem = problem()](const std::size_t eIdx)
            {
                const auto element = gridDiscretization.element(eIdx);
                const auto elemDisc = localView(gridDiscretization).bindElement(element);
                const auto elemSol = elementSolution(element, sol, gridDiscretization);

                std::vector<Variables> variables(Dumux::Detail::LocalDofs::numLocalDofs(elemDisc));
                for (const auto& localDof : localDofs(elemDisc))
                    variables[localDof.index()].update(elemSol, problem, elemDisc, ipData(elemDisc, localDof));

                (*ipDataCache)[eIdx].update(problem, element, elemDisc, variables);
            });
            ipDataCache_ = std::move(ipDataCache);
        }
    }

    template<class GridDiscretization, class SolutionVector>
    void update(const GridDiscretization& gridDiscretization, const SolutionVector& sol) {}

    const InterpolationPointData& scvfCache(std::size_t eIdx, std::size_t scvfIdx, std::size_t qpIdx) const
    { return (*ipDataCache_)[eIdx].scvfCache(scvfIdx, qpIdx); }

    const auto& cache(std::size_t eIdx) const
    { return (*ipDataCache_)[eIdx]; }

    const Problem& problem() const
    { return *problemPtr_;}

private:
    using ElementInterpolationPointDataCache = Dumux::Detail::CVFE::ScvfInterpolationPointDataCache<InterpolationPointData>;

    const Problem* problemPtr_;
    std::shared_ptr<const std::vector<ElementInterpolationPointDataCache>> ipDataCache_;
};

} // end namespace Dumux

#endif
