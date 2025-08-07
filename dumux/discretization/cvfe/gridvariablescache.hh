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
#ifndef DUMUX_DISCRETIZATION_CVFE_GRID_LOCAL_VARIABLES_HH
#define DUMUX_DISCRETIZATION_CVFE_GRID_LOCAL_VARIABLES_HH

#include <vector>
#include <type_traits>

#include <dumux/parallel/parallel_for.hh>

#include <dumux/common/concepts/localdofs_.hh>

// make the local view function available whenever we use this class
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/cvfe/elementvariables.hh>
#include <dumux/discretization/cvfe/elementsolution.hh>

namespace Dumux::Detail::CVFE {

template<class P, class V>
struct CVFEDefaultGridVariablesCacheTraits
{
    using Problem = P;
    using Variables = V;

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

    //! TODO: Delete after new variables concept has been implemented
    //! this is currently needed to support old interface
    using VolumeVariables = Variables;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = true;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    //! export the type of the mutable local view
    using MutableLocalView = LocalView::MutableView;

    CVFEGridVariablesCache(const Problem& problem) : problemPtr_(&problem) {}

    template<class GridGeometry, class SolutionVector>
    void update(const GridGeometry& gridGeometry, const SolutionVector& sol)
    {
        variables_.resize(gridGeometry.gridView().size(0));
        Dumux::parallelFor(gridGeometry.gridView().size(0), [&, &problem = problem()](const std::size_t eIdx)
        {
            const auto element = gridGeometry.element(eIdx);
            const auto fvGeometry = localView(gridGeometry).bindElement(element);

            // get the element solution
            auto elemSol = elementSolution(element, sol, gridGeometry);

            variables_[eIdx].resize(Dumux::Detail::LocalDofs::numLocalDofs(fvGeometry));
            for (const auto& localDof : localDofs(fvGeometry))
                variables_[eIdx][localDof.index()].update(elemSol, problem, fvGeometry, ipData(fvGeometry, localDof));
        });
    }

    //! TODO: Rename the functions below after new variables concept has been implemented
    //! this is currently needed to support old interface
    template<class ScvOrLocalDof, typename std::enable_if_t<!std::is_integral<ScvOrLocalDof>::value, int> = 0>
    const Variables& volVars(const ScvOrLocalDof& scvOrLocalDof) const
    {
        if constexpr (Concept::LocalDof<ScvOrLocalDof>)
            return variables_[scvOrLocalDof.elementIndex()][scvOrLocalDof.index()];
        else
            return variables_[scvOrLocalDof.elementIndex()][scvOrLocalDof.indexInElement()];
    }

    template<class ScvOrLocalDof, typename std::enable_if_t<!std::is_integral<ScvOrLocalDof>::value, int> = 0>
    Variables& volVars(const ScvOrLocalDof& scvOrLocalDof)
    {
        if constexpr (Concept::LocalDof<ScvOrLocalDof>)
            return variables_[scvOrLocalDof.elementIndex()][scvOrLocalDof.index()];
        else
            return variables_[scvOrLocalDof.elementIndex()][scvOrLocalDof.indexInElement()];
    }

    const Variables& volVars(const std::size_t eIdx, const std::size_t localIdx) const
    { return variables_[eIdx][localIdx]; }

    Variables& volVars(const std::size_t eIdx, const std::size_t localIdx)
    { return variables_[eIdx][localIdx]; }

    const Problem& problem() const
    { return *problemPtr_; }

private:
    const Problem* problemPtr_;
    std::vector<std::vector<Variables>> variables_;
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

    //! TODO: Delete after new variables concept has been implemented
    //! this is currently needed to support old interface
    using VolumeVariables = Variables;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = false;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    //! export the type of the mutable local view
    using MutableLocalView = LocalView::MutableView;

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
