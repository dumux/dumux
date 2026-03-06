// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FEDiscretization
 * \brief The grid local variables class for finite element methods
 */
#ifndef DUMUX_DISCRETIZATION_FE_GRID_LOCAL_VARIABLES_HH
#define DUMUX_DISCRETIZATION_FE_GRID_LOCAL_VARIABLES_HH

#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>
#include <ranges>
#include <memory>

#include <dumux/parallel/parallel_for.hh>

#include <dumux/common/concepts/localdofs_.hh>

// make the local view function available whenever we use this class
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/cvfe/elementsolution.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>
#include "elementvariables.hh"

namespace Dumux::Experimental::FE {

template<class P, class V, class IPD>
struct FEDefaultGridVariablesCacheTraits
{
    using Problem = P;
    using Variables = V;
    using InterpolationPointData = IPD;

    template<class GridVariablesCache, bool cachingEnabled>
    using LocalView = FEElementVariables<GridVariablesCache, cachingEnabled>;
};

/*!
 * \ingroup FEDiscretization
 * \brief Base class for the grid local variables
 */
template<class Traits, bool enableCaching>
class FEGridVariablesCache;

// specialization in case of storing the local variables
template<class Traits>
class FEGridVariablesCache<Traits, /*cachingEnabled*/true>
{
    using ThisType = FEGridVariablesCache<Traits, true>;
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

    FEGridVariablesCache(const Problem& problem) : problemPtr_(&problem) {}

    template<class GridGeometry, class SolutionVector>
    void init(const GridGeometry& gridGeometry, const SolutionVector& sol)
    {
        variables_.resize(gridGeometry.gridView().size(0));
        ipDataCache_ = std::make_shared<InterpolationPointDataCache>();
        ipDataCache_->resize(gridGeometry.gridView().size(0));

        Dumux::parallelFor(gridGeometry.gridView().size(0), [&, &problem = problem()](const std::size_t eIdx)
        {
            const auto element = gridGeometry.element(eIdx);
            const auto elemDisc = localView(gridGeometry).bindElement(element);

            // get the element solution
            auto elemSol = elementSolution(element, sol, gridGeometry);

            variables_[eIdx].resize(Dumux::Detail::LocalDofs::numLocalDofs(elemDisc));
            for (const auto& localDof : localDofs(elemDisc))
                variables_[eIdx][localDof.index()].update(elemSol, problem, elemDisc, ipData(elemDisc, localDof));

            ipDataCache_->update(problem, element, elemDisc, variables_[eIdx]);
        });
    }

    template<class GridGeometry, class SolutionVector>
    void update(const GridGeometry& gridGeometry, const SolutionVector& sol)
    {
        if constexpr (InterpolationPointData::isSolDependent)
        {
            auto newIpDataCache = std::make_shared<InterpolationPointDataCache>();
            newIpDataCache->resize(gridGeometry.gridView().size(0));

            Dumux::parallelFor(gridGeometry.gridView().size(0), [&, &problem = problem(), newIpDataCache](const std::size_t eIdx)
            {
                const auto element = gridGeometry.element(eIdx);
                const auto elemDisc = localView(gridGeometry).bindElement(element);

                // get the element solution
                auto elemSol = elementSolution(element, sol, gridGeometry);

                for (const auto& localDof : localDofs(elemDisc))
                    variables_[eIdx][localDof.index()].update(elemSol, problem, elemDisc, ipData(elemDisc, localDof));

                newIpDataCache->update(problem, element, elemDisc, variables_[eIdx]);
            });

            ipDataCache_ = std::move(newIpDataCache);
        }
        else
        {
            Dumux::parallelFor(gridGeometry.gridView().size(0), [&, &problem = problem()](const std::size_t eIdx)
            {
                const auto element = gridGeometry.element(eIdx);
                const auto elemDisc = localView(gridGeometry).bindElement(element);

                // get the element solution
                auto elemSol = elementSolution(element, sol, gridGeometry);

                for (const auto& localDof : localDofs(elemDisc))
                    variables_[eIdx][localDof.index()].update(elemSol, problem, elemDisc, ipData(elemDisc, localDof));
            });
        }
    }

    template<class LocalDof>
    const Variables& variables(const LocalDof& localDof) const
    { return variables_[localDof.elementIndex()][localDof.index()]; }

    template<class LocalDof>
    Variables& variables(const LocalDof& localDof)
    { return variables_[localDof.elementIndex()][localDof.index()]; }

    const Variables& variables(const std::size_t eIdx, const std::size_t localIdx) const
    { return variables_[eIdx][localIdx]; }

    Variables& variables(const std::size_t eIdx, const std::size_t localIdx)
    { return variables_[eIdx][localIdx]; }

    const InterpolationPointData& elementCache(std::size_t eIdx, std::size_t qpIdx) const
    { return ipDataCache_->elementCache(eIdx, qpIdx); }

    InterpolationPointData& elementCache(std::size_t eIdx, std::size_t qpIdx)
    { return ipDataCache_->elementCache(eIdx, qpIdx); }

    const InterpolationPointData& boundaryIntersectionCache(std::size_t eIdx, int intersectionIdx, std::size_t qpIdx) const
    { return ipDataCache_->boundaryIntersectionCache(eIdx, intersectionIdx, qpIdx); }

    InterpolationPointData& boundaryIntersectionCache(std::size_t eIdx, int intersectionIdx, std::size_t qpIdx)
    { return ipDataCache_->boundaryIntersectionCache(eIdx, intersectionIdx, qpIdx); }

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
            elementCache_.resize(numElements);
            boundaryIntersectionCache_.resize(numElements);
        }

        template<class Problem, class ElementDiscretization, class ElementVariables>
        void update(const Problem& problem,
                    const typename ElementDiscretization::Element& element,
                    const ElementDiscretization& elemDisc,
                    const ElementVariables& elemVars)
        {
            const auto eIdx = elemDisc.gridGeometry().elementMapper().index(element);
            updateElementCache_(problem, eIdx, element, elemDisc, elemVars);
        }

        // access operator
        const InterpolationPointData& elementCache(std::size_t eIdx, std::size_t qpIdx) const
        { return elementCache_[eIdx][qpIdx]; }

        // access operator
        InterpolationPointData& elementCache(std::size_t eIdx, std::size_t qpIdx)
        { return elementCache_[eIdx][qpIdx]; }

        // access operator
        const InterpolationPointData& boundaryIntersectionCache(std::size_t eIdx, int intersectionIdx, std::size_t qpIdx) const
        { return (*boundaryIntersectionCache_[eIdx]).at(intersectionIdx)[qpIdx]; }

        // access operator
        InterpolationPointData& boundaryIntersectionCache(std::size_t eIdx, int intersectionIdx, std::size_t qpIdx)
        { return (*boundaryIntersectionCache_[eIdx])[intersectionIdx][qpIdx]; }

    private:
        template<class Problem, class ElementDiscretization, class ElementVariables>
        void updateElementCache_(const Problem& problem,
                                 const std::size_t eIdx,
                                 const typename ElementDiscretization::Element& element,
                                 const ElementDiscretization& elemDisc,
                                 const ElementVariables& elemVars)
        {
            const auto elemQuadRule = Dumux::CVFE::quadratureRule(elemDisc, element);
            elementCache_[eIdx].resize(std::ranges::size(elemQuadRule));
            for (const auto& qpData : elemQuadRule)
                elementCache_[eIdx][qpData.ipData().qpIndex()].update(problem, element, elemDisc, elemVars, qpData.ipData());

            if (!boundaryIntersectionCache_[eIdx])
                boundaryIntersectionCache_[eIdx] = std::make_unique<std::unordered_map<int, std::vector<InterpolationPointData>>>();
            else
                boundaryIntersectionCache_[eIdx]->clear();

            for (const auto& intersection : intersections(elemDisc.gridGeometry().gridView(), element))
            {
                if (intersection.boundary())
                {
                    const auto intersectionIndex = intersection.indexInInside();
                    auto& boundaryCache = (*boundaryIntersectionCache_[eIdx])[intersectionIndex];
                    const auto quadRule = Dumux::CVFE::quadratureRule(elemDisc, intersection);
                    boundaryCache.resize(std::ranges::size(quadRule));
                    for (const auto& qpData : quadRule)
                        boundaryCache[qpData.ipData().qpIndex()].update(problem,
                                                                        element,
                                                                        elemDisc,
                                                                        elemVars,
                                                                        qpData.ipData());
                }
            }
        }

        std::vector<std::vector<InterpolationPointData>> elementCache_; //! storage per element quadrature points
        std::vector<std::unique_ptr<std::unordered_map<int, std::vector<InterpolationPointData>>>> boundaryIntersectionCache_; //! storage per element/boundary intersection/qp
    };

    const Problem* problemPtr_;
    std::vector<std::vector<Variables>> variables_;
    std::shared_ptr<InterpolationPointDataCache> ipDataCache_;
};

// Specialization when the current local variables are not stored
template<class Traits>
class FEGridVariablesCache<Traits, /*cachingEnabled*/false>
{
    using ThisType = FEGridVariablesCache<Traits, false>;

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

    FEGridVariablesCache(const Problem& problem) : problemPtr_(&problem) {}

    template<class GridGeometry, class SolutionVector>
    void update(const GridGeometry& gridGeometry, const SolutionVector& sol) {}

    const Problem& problem() const
    { return *problemPtr_;}

private:
    const Problem* problemPtr_;
};

} // end namespace Dumux::Experimental::CVFE

#endif
