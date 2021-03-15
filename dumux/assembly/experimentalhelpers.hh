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
#include <dumux/discretization/concepts.hh>

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

// get the Dirichlet values for a boundary entity from a user problem
template<class Problem, class Element, class BoundaryEntity, class ElemSol>
decltype(auto) getDirichletValues(const Problem& problem,
                                  const Element& element,
                                  const BoundaryEntity& boundaryEntity,
                                  const ElemSol& elemSol)
{
    // given element solution models the experimental elemsol concept
    if constexpr (Dune::models<Concept::ElementSolution, ElemSol>())
    {
        if constexpr (ProblemTraits<Problem>::hasTransientDirichletInterface)
            return problem.dirichlet(element, boundaryEntity, elemSol.timeLevel());
        else
            return problem.dirichlet(element, boundaryEntity);
    }

    // current default style (elemSolState = element solution, not passed to user interface)
    else
        return problem.dirichlet(element, boundaryEntity);
}

} // end namespace CompatibilityHelpers
} // end namespace Dumux::Experimental

#endif
