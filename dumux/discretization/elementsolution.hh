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
 * \brief Element solution classes and factory functions
 */
#ifndef DUMUX_DISCRETIZATION_ELEMENT_SOLUTION_HH
#define DUMUX_DISCRETIZATION_ELEMENT_SOLUTION_HH

#include <type_traits>

#include <dune/common/concept.hh>
#include <dumux/discretization/concepts.hh>

#include <dumux/discretization/cellcentered/elementsolution.hh>
#include <dumux/discretization/box/elementsolution.hh>
#include <dumux/discretization/staggered/elementsolution.hh>

namespace Dumux {

struct EmptyElementSolution {};

namespace Experimental {

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

// hide deduction guides from doxygen
#ifndef DOXYGEN
template<class B, class TL> ElementSolution(const B& b, const TL& t) -> ElementSolution<B, TL>;
template<class B, class TL> ElementSolution(B&& b, const TL& t) -> ElementSolution<B, TL>;
#endif // DOXYGEN

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
    return ElementSolution{std::move(elemSol), solState.timeLevel()};
}

} // end namespace Experimental
} // end namespace Dumux

#endif
