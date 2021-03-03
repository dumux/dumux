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
 * \ingroup Discretization
 * \brief TODO: Doc me
 */
#ifndef DUMUX_SOLUTION_STATE_HH
#define DUMUX_SOLUTION_STATE_HH

#include <dune/common/concept.hh>
#include <dumux/discretization/elementsolution.hh>

namespace Dumux {
namespace Experimental {

namespace Concept {

//! Concept for solution states
struct SolutionState
{
  template<class T>
  auto require(const T& t) -> decltype(
    t.dofs(),
    t.timeLevel()
  );
};

//! Concept for element solution states
struct ElementSolutionState
{
  template<class T>
  auto require(const T& t) -> decltype(
    t.elementSolution(),
    t.timeLevel()
  );
};

} // end namespace Concept

/*!
 * \ingroup Discretization
 * \brief TODO: Doc me
 */
template<class SV, class TL>
class SolutionState
{
public:
    using SolutionVector = SV;
    using TimeLevel = TL;

    //! Constructor
    SolutionState(const SolutionVector* dofs,
                  const TimeLevel* timeLevel)
    : dofs_(dofs)
    , timeLevel_(timeLevel)
    {}

    //! return the solution vector
    const SolutionVector& dofs() const
    { return *dofs_; }

    //! return the time level
    const TimeLevel& timeLevel() const
    { return *timeLevel_; }

private:
    const SolutionVector* dofs_;
    const TimeLevel* timeLevel_;
};

/*!
 * \ingroup Discretization
 * \brief TODO: Doc me
 */
template<class ES, class TL>
class ElementSolutionState
{
public:
    using ElementSolution = ES;
    using TimeLevel = TL;

    //! Constructor
    ElementSolutionState(ElementSolution&& es,
                         const TimeLevel* timeLevel)
    : elemSol_(std::move(es))
    , timeLevel_(timeLevel)
    {}

    //! return the element solution vector
    const ElementSolution& elementSolution() const
    { return elemSol_; }

    //! return the time level
    const TimeLevel& timeLevel() const
    { return *timeLevel_; }

private:
    ElementSolution elemSol_;
    const TimeLevel* timeLevel_;
};

/*!
 * \ingroup Discretization
 * \brief TODO: Doc me
 */
template<class Element, class SolutionState, class GridGeometry>
auto makeElementSolutionState(const Element& element,
                              const SolutionState& state,
                              const GridGeometry& gridGeom)
{
    static_assert(Dune::models<Concept::SolutionState, SolutionState>(),
                  "Provided type does not fulfill solution state interface");

    auto elemSol = elementSolution(element, state.dofs(), gridGeom);
    return ElementSolutionState{std::move(elemSol), &state.timeLevel()};
}

} // end namespace Experimental
} // end namespace Dumux

#endif
