// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CVFEDiscretization
 * \brief The grid local variables class for control-volume finite element methods
 */
#ifndef DUMUX_DISCRETIZATION_CVFE_GRID_LOCAL_VARIABLES__HH
#define DUMUX_DISCRETIZATION_CVFE_GRID_LOCAL_VARIABLES__HH

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
 * \brief Base class for the grid local variables
 */
template<class Traits, bool enableCaching>
class CVFEGridVariablesCache;

// specialization in case of storing the local variables
template<class Traits>
class CVFEGridVariablesCache<Traits, /*cachingEnabled*/true>
{
    using ThisType = CVFEGridVariablesCache<Traits, true>;
public:
    //! export the problem type
    using Problem = typename Traits::Problem;

    //! export the variables type
    using Variables = typename Traits::Variables;

    //! export interpolation point data type
    using InterpolationPointData = typename Traits::InterpolationPointData;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = true;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    //! export the type of the mutable local view
    using MutableLocalView = LocalView::MutableView;

    CVFEGridVariablesCache(const Problem& problem) : problemPtr_(&problem) {}

    template<class GridGeometry, class SolutionVector>
    void init(const GridGeometry& gridGeometry, const SolutionVector& sol)
    {
        variables_.resize(gridGeometry.gridView().size(0));
        ipDataCache_ = std::make_shared<InterpolationPointDataCache>();
        ipDataCache_->resize(gridGeometry.gridView().size(0));

        Dumux::parallelFor(gridGeometry.gridView().size(0), [&, &problem = problem()](const std::size_t eIdx)
        {
            const auto element = gridGeometry.element(eIdx);
            const auto fvGeometry = localView(gridGeometry).bindElement(element);

            // get the element solution
            auto elemSol = elementSolution(element, sol, gridGeometry);

            variables_[eIdx].resize(Dumux::Detail::LocalDofs::numLocalDofs(fvGeometry));
            for (const auto& localDof : localDofs(fvGeometry))
                variables_[eIdx][localDof.index()].update(elemSol, problem, fvGeometry, ipData(fvGeometry, localDof));

            ipDataCache_->update(problem, element, fvGeometry, variables_[eIdx]);
        });
    }

    template<class GridGeometry, class SolutionVector>
    void update(const GridGeometry& gridGeometry, const SolutionVector& sol)
    {
        if constexpr (InterpolationPointData::isSolDependent)
        {
            auto newIpDataCache = std::make_shared<InterpolationPointDataCache>(*ipDataCache_);

            Dumux::parallelFor(gridGeometry.gridView().size(0), [&, &problem = problem(), newIpDataCache](const std::size_t eIdx)
            {
                const auto element = gridGeometry.element(eIdx);
                const auto fvGeometry = localView(gridGeometry).bindElement(element);

                // get the element solution
                auto elemSol = elementSolution(element, sol, gridGeometry);

                for (const auto& localDof : localDofs(fvGeometry))
                    variables_[eIdx][localDof.index()].update(elemSol, problem, fvGeometry, ipData(fvGeometry, localDof));

                newIpDataCache->update(problem, element, fvGeometry, variables_[eIdx]);
            });

            ipDataCache_ = std::move(newIpDataCache);
        }
        else
        {
            Dumux::parallelFor(gridGeometry.gridView().size(0), [&, &problem = problem()](const std::size_t eIdx)
            {
                const auto element = gridGeometry.element(eIdx);
                const auto fvGeometry = localView(gridGeometry).bindElement(element);

                // get the element solution
                auto elemSol = elementSolution(element, sol, gridGeometry);

                for (const auto& localDof : localDofs(fvGeometry))
                    variables_[eIdx][localDof.index()].update(elemSol, problem, fvGeometry, ipData(fvGeometry, localDof));
            });
        }
    }

    template<Dumux::Concept::LocalDof LocalDof>
    const Variables& variables(const LocalDof& localDof) const
    { return variables_[localDof.elementIndex()][localDof.index()]; }

    template<Dumux::Concept::LocalDof LocalDof>
    Variables& variables(const LocalDof& localDof)
    { return variables_[localDof.elementIndex()][localDof.index()]; }

    const Variables& variables(const std::size_t eIdx, const std::size_t localIdx) const
    { return variables_[eIdx][localIdx]; }

    Variables& variables(const std::size_t eIdx, const std::size_t localIdx)
    { return variables_[eIdx][localIdx]; }

    const InterpolationPointData& cache(std::size_t eIdx, std::size_t scvfIdx, std::size_t qpIdx) const
    { return ipDataCache_->cache(eIdx, scvfIdx, qpIdx); }

    InterpolationPointData& cache(std::size_t eIdx, std::size_t scvfIdx, std::size_t qpIdx)
    { return ipDataCache_->cache(eIdx, scvfIdx, qpIdx); }

    const Problem& problem() const
    { return *problemPtr_; }

private:
    class InterpolationPointDataCache
    {
    public:
        InterpolationPointDataCache()
        {}

        void resize(const std::size_t numElements)
        {
            scvfCache_.resize(numElements);
            qpsOffset_.resize(numElements);
        }

        template<class Problem, class FVElementGeometry, class ElementVariables>
        void update(const Problem& problem,
                    const typename FVElementGeometry::Element& element,
                    const FVElementGeometry& fvGeometry,
                    const ElementVariables& elemVars)
        {
            const auto eIdx = fvGeometry.gridGeometry().elementMapper().index(element);
            updateElementCache_(problem, eIdx, element, fvGeometry, elemVars);
        }

        // access operator
        const InterpolationPointData& cache(std::size_t eIdx, std::size_t scvfIdx, std::size_t qpIdx) const
        { return scvfCache_[eIdx][qpsOffset_[eIdx][scvfIdx] + qpIdx]; }

        // access operator
        InterpolationPointData& cache(std::size_t eIdx, std::size_t scvfIdx, std::size_t qpIdx)
        { return scvfCache_[eIdx][qpsOffset_[eIdx][scvfIdx] + qpIdx]; }

    private:
        template<class Problem, class FVElementGeometry, class ElementVariables>
        void updateElementCache_(const Problem& problem,
                                 const std::size_t eIdx,
                                 const typename FVElementGeometry::Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElementVariables& elemVars)
        {
            qpsOffset_[eIdx].resize(fvGeometry.numScvf() + 1, 0);
            for (const auto& scvf : scvfs(fvGeometry))
            {
                const auto numQps = std::ranges::size(Dumux::CVFE::quadratureRule(fvGeometry, scvf));
                qpsOffset_[eIdx][scvf.index() + 1] = numQps;
            }
            for (std::size_t i = 2; i < qpsOffset_[eIdx].size(); ++i)
                qpsOffset_[eIdx][i] += qpsOffset_[eIdx][i-1];

            scvfCache_[eIdx].resize(qpsOffset_[eIdx].back());
            for (const auto& scvf : scvfs(fvGeometry))
            {
                for (const auto& qpData : Dumux::CVFE::quadratureRule(fvGeometry, scvf))
                    cache(eIdx, qpData.ipData().scvfIndex(), qpData.ipData().qpIndex()).update(problem,
                                                                                                element,
                                                                                                fvGeometry,
                                                                                                elemVars,
                                                                                                qpData.ipData());
            }
        }

        std::vector<std::vector<InterpolationPointData>> scvfCache_; //! storage per element
        std::vector<std::vector<std::size_t>> qpsOffset_; //! offset[eIdx][scvfIdx] -> start index in scvfCache_[eIdx]
    };

    const Problem* problemPtr_;
    std::vector<std::vector<Variables>> variables_;
    std::shared_ptr<InterpolationPointDataCache> ipDataCache_;
};

// Specialization when the current local variables are not stored
template<class Traits>
class CVFEGridVariablesCache<Traits, /*cachingEnabled*/false>
{
    using ThisType = CVFEGridVariablesCache<Traits, false>;

public:
    //! export the problem type
    using Problem = typename Traits::Problem;

    //! export the variables type
    using Variables = typename Traits::Variables;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = false;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    //! export the type of the mutable local view
    using MutableLocalView = LocalView::MutableView;

    //! export interpolation point data type
    using InterpolationPointData = typename Traits::InterpolationPointData;

    CVFEGridVariablesCache(const Problem& problem) : problemPtr_(&problem) {}

    template<class GridGeometry, class SolutionVector>
    void update(const GridGeometry& gridGeometry, const SolutionVector& sol) {}

    const Problem& problem() const
    { return *problemPtr_;}

private:
    const Problem* problemPtr_;
};

} // end namespace Dumux

#endif
