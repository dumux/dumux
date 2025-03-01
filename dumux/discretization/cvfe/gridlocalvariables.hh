// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
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

// make the local view function available whenever we use this class
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/cvfe/elementlocalvariables.hh>
#include <dumux/discretization/cvfe/elementsolution.hh>

namespace Dumux {

template<class P, class VV>
struct CVFEDefaultGridLocalVariablesTraits
{
    using Problem = P;
    using LocalVariables = VV;

    template<class GridLocalVariables, bool cachingEnabled>
    using LocalView = CVFEElementLocalVariables<GridLocalVariables, cachingEnabled>;
};

/*!
 * \ingroup CVFEDiscretization
 * \brief Base class for the grid local variables
 */
template<class Traits, bool enableCaching>
class CVFEGridLocalVariables;

// specialization in case of storing the local variables
template<class Traits>
class CVFEGridLocalVariables<Traits, /*cachingEnabled*/true>
{
    using ThisType = CVFEGridLocalVariables<Traits, true>;

public:
    //! export the problem type
    using Problem = typename Traits::Problem;

    //! export the local variables type
    using LocalVariables = typename Traits::LocalVariables;

    //! export the volume variables type
    using VolumeVariables [[deprecated("Use LocalVariables instead. Will be removed after release 3.10.")]] = LocalVariables;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = true;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    CVFEGridLocalVariables(const Problem& problem) : problemPtr_(&problem) {}

    template<class GridGeometry, class SolutionVector>
    void update(const GridGeometry& gridGeometry, const SolutionVector& sol)
    {
        localVariables_.resize(gridGeometry.gridView().size(0));
        Dumux::parallelFor(gridGeometry.gridView().size(0), [&, &problem = problem()](const std::size_t eIdx)
        {
            const auto element = gridGeometry.element(eIdx);
            const auto fvGeometry = localView(gridGeometry).bindElement(element);

            // get the element solution
            auto elemSol = elementSolution(element, sol, gridGeometry);

            localVariables_[eIdx].resize(Detail::numLocalDofs(fvGeometry));
            for (const auto& localDof : localDofs(fvGeometry))
                localVariables_[eIdx][localDof.index()].update(elemSol, problem, fvGeometry, localDof);
        });
    }

    template<class ScvOrLocalDof, typename std::enable_if_t<!std::is_integral<ScvOrLocalDof>::value, int> = 0>
    const LocalVariables& localVars(const ScvOrLocalDof& scvOrLocalDof) const
    {
        if constexpr (Detail::isLocalDofType<ScvOrLocalDof>())
            return localVariables_[scvOrLocalDof.elementIndex()][scvOrLocalDof.index()];
        else
            return localVariables_[scvOrLocalDof.elementIndex()][scvOrLocalDof.indexInElement()];
    }

    template<class ScvOrLocalDof, typename std::enable_if_t<!std::is_integral<ScvOrLocalDof>::value, int> = 0>
    LocalVariables& localVars(const ScvOrLocalDof& scvOrLocalDof)
    {
        if constexpr (Detail::isLocalDofType<ScvOrLocalDof>())
            return localVariables_[scvOrLocalDof.elementIndex()][scvOrLocalDof.index()];
        else
            return localVariables_[scvOrLocalDof.elementIndex()][scvOrLocalDof.indexInElement()];
    }

    const LocalVariables& localVars(const std::size_t eIdx, const std::size_t localIdx) const
    { return localVariables_[eIdx][localIdx]; }

    LocalVariables& localVars(const std::size_t eIdx, const std::size_t localIdx)
    { return localVariables_[eIdx][localIdx]; }

    const Problem& problem() const
    { return *problemPtr_; }

private:
    const Problem* problemPtr_;
    std::vector<std::vector<LocalVariables>> localVariables_;
};

// Specialization when the current local variables are not stored
template<class Traits>
class CVFEGridLocalVariables<Traits, /*cachingEnabled*/false>
{
    using ThisType = CVFEGridLocalVariables<Traits, false>;

public:
    //! export the problem type
    using Problem = typename Traits::Problem;

    //! export the local variables type
    using LocalVariables = typename Traits::LocalVariables;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = false;

    //! export the volume variables type
    using VolumeVariables [[deprecated("Use LocalVariables instead. Will be removed after release 3.10.")]] = LocalVariables;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    CVFEGridLocalVariables(const Problem& problem) : problemPtr_(&problem) {}

    template<class GridGeometry, class SolutionVector>
    void update(const GridGeometry& gridGeometry, const SolutionVector& sol) {}

    const Problem& problem() const
    { return *problemPtr_;}

private:
    const Problem* problemPtr_;
};

} // end namespace Dumux

#endif
