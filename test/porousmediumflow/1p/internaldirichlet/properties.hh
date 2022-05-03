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
 * \ingroup OnePTests
 * \brief A test for internal Dirichlet constraints
 */

#ifndef DUMUX_INCOMPRESSIBLE_ONEP_TEST_PROBLEM_INTERNAL_DIRICHLET_PROPERTIES_HH
#define DUMUX_INCOMPRESSIBLE_ONEP_TEST_PROBLEM_INTERNAL_DIRICHLET_PROPERTIES_HH

#include <test/porousmediumflow/1p/incompressible/properties.hh>
#include "problem.hh"

namespace Dumux::Properties {
// Create new type tags
namespace TTag {
template<class T> struct OnePInternalDirichlet {};
template<class T> struct OnePInternalDirichletTpfa { using InheritsFrom = std::tuple<OnePInternalDirichlet<T>, OnePIncompressibleTpfa<T>>; using Traits = T; };
template<class T> struct OnePInternalDirichletBox { using InheritsFrom = std::tuple<OnePInternalDirichlet<T>, OnePIncompressibleBox<T>>; using Traits = T; };
} // end namespace TTag

// Set the problem type
template<class TypeTag, class T>
struct Problem<TypeTag, TTag::OnePInternalDirichlet<T>>
{ using type = OnePTestProblemInternalDirichlet<TypeTag>; };

} // end namespace Dumux::Properties

#endif
