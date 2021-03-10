// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup Assembly
 * \brief Helper classes and functions to be used in places where compatibility
 *        between the current default and the new experimental assembly is seeked.
 */
#ifndef DUMUX_ASSEMBLY_EXPERIMENTAL_HELPERS_HH
#define DUMUX_ASSEMBLY_EXPERIMENTAL_HELPERS_HH

#include <vector>
#include <utility>
#include <type_traits>

#include <dune/common/fvector.hh>

#include <dumux/common/typetraits/problem.hh>
#include <dumux/discretization/solutionstate.hh>
#include <dumux/discretization/fvgridvariables.hh>
#include <dumux/discretization/elementsolution.hh>

namespace Dumux::Experimental {
namespace CompatibilityHelpers {

// implements necessary interfaces for compatibility
template<class ElemVolVars, class ElemFluxVarsCache>
class GridVariablesFacade
{
    using P = typename ElemVolVars::GridVolumeVariables::Problem;
public:
    using GridVolumeVariables = typename ElemVolVars::GridVolumeVariables;
    using GridFluxVariablesCache = typename ElemFluxVarsCache::GridFluxVariablesCache;
    using GridGeometry = typename ProblemTraits<P>::GridGeometry;
};

// implements necessary interfaces for compatibility
template<class ElemVolVars, class ElemFluxVarsCache>
class ElemVariablesFacade
{
    const ElemVolVars* elemVolVars_;
    const ElemFluxVarsCache* elemFluxVarsCache_;

public:

    using GridVariables = GridVariablesFacade<ElemVolVars, ElemFluxVarsCache>;

    ElemVariablesFacade(const ElemVolVars* elemVolVars,
                        const ElemFluxVarsCache* elemFluxVarsCache)
    : elemVolVars_(elemVolVars)
    , elemFluxVarsCache_(elemFluxVarsCache)
    {}

    const ElemVolVars& elemVolVars() const { return *elemVolVars_; }
    const ElemFluxVarsCache& elemFluxVarsCache() const { return *elemFluxVarsCache_; }
};

// helper function to construct an element variables facade from given local views
template<class ElemVolVars, class ElemFluxVarsCache>
ElemVariablesFacade<ElemVolVars, ElemFluxVarsCache>
makeElemVariablesFacade(const ElemVolVars& elemVolVars,
                        const ElemFluxVarsCache& elemFluxVarsCache)
{ return {&elemVolVars, &elemFluxVarsCache}; }

// Get an element solution state from a global state.
// If a solution vector is passed (current default style),
// then an element solution is returned, whereas if an actual
// solution state is given (new concept), we return an elem sol state.
template<class Element, class SolState, class GridGeometry>
auto getElemSolState(const Element& element,
                     const SolState& solState,
                     const GridGeometry& gridGeometry)
{
    if constexpr (Dune::models<Concept::SolutionState, SolState>())
        return makeElementSolutionState(element, solState, gridGeometry);
    else
        return elementSolution(element, solState, gridGeometry);
}

// Get an element solution state from grid variables.
// If grid variables with the experimental layout are passed,
// we return an actual element solution state. Otherwise,
// we return an empty element solution.
template<class Element, class GridVariables>
auto getElemSolState(const Element& element,
                     const GridVariables& gridVariables)
{
    if constexpr (areExperimentalGridVars<GridVariables>)
        return makeElementSolutionState(element, gridVariables);
    else
    {
        using DummySolVec = std::vector<Dune::FieldVector<double, 1>>;
        using ElemSol = decltype(elementSolution(element,
                                                 std::declval<DummySolVec>(),
                                                 gridVariables.gridGeometry()));
        return ElemSol{};
    }
}

// get the Dirichlet values for a boundary entity from a user problem
template<class Problem, class Element, class BoundaryEntity, class ElemSolState>
decltype(auto) getDirichletValues(const Problem& problem,
                                  const Element& element,
                                  const BoundaryEntity& boundaryEntity,
                                  const ElemSolState& elemSolState)
{
    // given element solution state models an actual state
    if constexpr (Dune::models<Concept::ElementSolutionState, ElemSolState>())
        return problem.dirichlet(element, boundaryEntity, elemSolState);

    // current default style (elemSolState = element solution, not passed to user interface)
    else
        return problem.dirichlet(element, boundaryEntity);
}

} // end namespace CompatibilityHelpers
} // end namespace Dumux::Experimental

#endif
