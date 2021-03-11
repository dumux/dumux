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
 * \brief State classes to represent the global and element-local states of the
 *        primary variables of a numerical model, consisting of the spatial degrees
 *        of freedom and potentially a corresponding time level.
 */
#ifndef DUMUX_SOLUTION_STATE_HH
#define DUMUX_SOLUTION_STATE_HH

#include <utility>
#include <type_traits>

#include <dune/common/concept.hh>
#include <dumux/discretization/elementsolution.hh>

namespace Dumux::Experimental {

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

//! Concept for element solutions
struct ElementSolution
{
    template<class T>
    auto require(const T& t) -> decltype(
        t.size(),
        t[0],
        t.timeLevel()
    );
};

} // end namespace Concept

/*!
 * \ingroup Discretization
 * \brief View on the state of a numerical solution, consisting of the spatial
 *        degrees of freedom and a corresponding time level.
 */
template<class SV, class TL>
class SolutionStateView
{
public:
    using SolutionVector = SV;
    using TimeLevel = TL;

    //! Constructor
    template<class SolutionVector, class TimeLevel>
    SolutionStateView(SolutionVector&& dofs,
                      TimeLevel&& timeLevel)
    : dofs_(dofs)
    , timeLevel_(timeLevel)
    {
        static_assert(std::is_lvalue_reference_v<SolutionVector>,
                      "Solution vector must be an lvalue reference");
        static_assert(std::is_lvalue_reference_v<TimeLevel>,
                      "Time level must be an lvalue reference");
    }

    //! Constructor from a solution state
    template<class SolState>
    SolutionStateView(const SolState& other)
    : SolutionStateView(other.dofs(), other.timeLevel())
    {
        static_assert(Dune::models<Concept::SolutionState, SolState>(),
                      "Given type does not model SolutionState concept");
    }

    //! return the solution vector
    const SolutionVector& dofs() const
    { return dofs_; }

    //! return the time level
    const TimeLevel& timeLevel() const
    { return timeLevel_; }

private:
    const SolutionVector& dofs_;
    const TimeLevel& timeLevel_;
};

/*!
 * \ingroup Discretization
 * \brief Make a solution state view from a solution state
 */
template<class SolutionState>
auto solutionStateView(const SolutionState& solState)
{
    static_assert(Dune::models<Concept::SolutionState, SolutionState>(),
                  "Given type does not model SolutionState concept");
    using SV = std::decay_t<decltype(solState.dofs())>;
    using TL = std::decay_t<decltype(solState.timeLevel())>;
    return SolutionStateView<SV, TL>(solState);
}

/*!
 * \ingroup Discretization
 * \brief State class to represent the element-local state of the primary
 *        variables of a numerical model, consisting of the spatial degrees
 *        of freedom and potentially a corresponding time level.
 * \note Here, we extend the concept of elementSolution to additionally carry
 *       time information. Therefore, for now we inherit from the given element
 *       solution, extending its interface
 */
template<class BaseElemSol, class TL>
class ElementSolution : public BaseElemSol
{
public:
    using TimeLevel = TL;

    //! Constructor
    template<class TheTimeLevel,
             std::enable_if_t<std::is_same_v<std::decay_t<TheTimeLevel>, TL>, int> = 0>
    ElementSolution(BaseElemSol&& base,
                    TheTimeLevel&& timeLevel)
    : BaseElemSol(std::forward<BaseElemSol>(base))
    , timeLevel_(timeLevel)
    {
        static_assert(std::is_lvalue_reference_v<TheTimeLevel>,
                      "Time level must be an lvalue reference");
    }

    //! return the time level
    const TimeLevel& timeLevel() const
    { return timeLevel_; }

private:
    const TimeLevel& timeLevel_;
};

/*!
 * \ingroup Discretization
 * \brief Make an element solution from a global solution state.
 */
template<class Element, class SolState, class GridGeometry,
         std::enable_if_t<Dune::models<Concept::SolutionState, SolState>(), int> = 0>
auto elementSolution(const Element& element,
                     const SolState& solState,
                     const GridGeometry& gridGeometry)
{
    auto elemSol = elementSolution(element, solState.dofs(), gridGeometry);

    using Base = std::decay_t<decltype(elemSol)>;
    using TL = std::decay_t<decltype(solState.timeLevel())>;
    return ElementSolution<Base, TL>{std::move(elemSol), solState.timeLevel()};
}

} // end namespace Dumux::Experimental

#endif
