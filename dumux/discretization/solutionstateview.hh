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
#ifndef DUMUX_DISCRETIZATION_SOLUTION_STATE_HH
#define DUMUX_DISCRETIZATION_SOLUTION_STATE_HH

#include <utility>
#include <type_traits>

#include <dune/common/concept.hh>
#include <dumux/discretization/concepts.hh>

namespace Dumux::Experimental {

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

// hide deduction guides from doxygen
#ifndef DOXYGEN
template<class SolState> SolutionStateView(const SolState& b)
    -> SolutionStateView<typename SolState::SolutionVector,
                         typename SolState::TimeLevel>;
#endif // DOXYGEN

} // end namespace Dumux::Experimental

#endif
