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
 * \brief Concepts related to discretization methods.
 */
#ifndef DUMUX_DISCRETIZATION_CONCEPT_HH
#define DUMUX_DISCRETIZATION_CONCEPT_HH

#include <dune/common/concept.hh>

namespace Dumux{

// experimental concepts
namespace Experimental::Concept {

//! Concept for states of discrete numerical solutions
struct SolutionState
{
    template<class T>
    auto require(const T& t) -> decltype(
        Dune::Concept::requireType<typename T::SolutionVector>(),
        Dune::Concept::requireType<typename T::TimeLevel>(),
        Dune::Concept::requireConvertible<typename T::SolutionVector>(t.dofs()),
        Dune::Concept::requireConvertible<typename T::TimeLevel>(t.timeLevel())
    );
};

//! Concept for element-local discrete numerical solutions
struct ElementSolution
{
    template<class T>
    auto require(const T& t) -> decltype(
        Dune::Concept::requireType<typename T::PrimaryVariables>(),
        Dune::Concept::requireType<typename T::TimeLevel>(),
        Dune::Concept::requireConvertible<std::size_t>(t.size()),
        Dune::Concept::requireConvertible<typename T::PrimaryVariables>(t[0]),
        Dune::Concept::requireConvertible<typename T::TimeLevel>(t.timeLevel())
    );
};

} // end namespace Experimental::Concept
} // end namespace Dumux

#endif
