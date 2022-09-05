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
 * \ingroup Common
 * \brief Interoperability functions for types that we support as
 *        variables of systems of equations.
 */
#ifndef DUMUX_COMMON_VARIABLES_HH
#define DUMUX_COMMON_VARIABLES_HH

#include <concepts>
#include <type_traits>
#include <utility>
#include <memory>

// include the default implementation
#include <dumux/experimental/new_assembly/dumux/common/defaultvariables.hh>

#include <dumux/experimental/new_assembly/dumux/common/concepts.hh>
#include <dumux/experimental/new_assembly/dumux/common/timelevel.hh>
#include <dumux/experimental/new_assembly/dumux/linear/system.hh>

namespace Dumux {
namespace Concepts {

//! Variables that expose dofs and allow updating them.
template<typename T>
concept VariablesWithDofs = requires(T& vars, const T& constVars) {
    typename T::Dofs;
    { constVars.dofs() } -> LinearSystem::Interoperable;
    { vars.update(constVars.dofs()) };
};

//! Variables can be either linear system types (vectors) or exposing dofs
template<typename T>
concept Variables = LinearSystem::Interoperable<T> or VariablesWithDofs<T>;

//! Concept for time-dependent variables
template<typename T>
concept TimeDependentVariables = VariablesWithDofs<T> and requires(T& vars, const T& constVars) {
        typename T::TimeLevel;
        { constVars.timeLevel() } -> std::convertible_to<typename T::TimeLevel>;
        { vars.update(constVars.dofs(), constVars.timeLevel()) };
        { vars.updateTime(constVars.timeLevel()) };
};

} // namespace Concepts


namespace Variables {

//! Get the degrees of freedom of variables
template<Concepts::VariablesWithDofs Variables>
decltype(auto) getDofs(const Variables& vars)
{ return vars.dofs(); }


//! Update variables with the given dofs
template<Concepts::VariablesWithDofs Variables>
void update(Variables& vars, const typename Variables::Dofs& dofs)
{ vars.update(dofs); }


//! Overload for linear system types
template<LinearSystem::Interoperable Dofs>
const Dofs& getDofs(const Dofs& dofs)
{ return dofs; }


//! Overload for linear system types
template<LinearSystem::Interoperable Dofs, typename NewDofs> requires(
    std::assignable_from<Dofs&, const NewDofs&>)
void update(Dofs& dofs, const NewDofs& newDofs)
{ dofs = newDofs; }


#ifndef DOXYGEN
namespace Detail {

template<typename T>
struct Dofs;

template<Concepts::VariablesWithDofs T>
struct Dofs<T> : public std::type_identity<typename T::Dofs> {};

template<LinearSystem::Interoperable T>
struct Dofs<T> : public std::type_identity<T> {};

} // namespace Detail
#endif // DOXYGEN


//! Convenience alias of the type of dofs used by variables
template<Concepts::Variables V>
using DofsType = typename Detail::Dofs<V>::type;


//! Convenience alias of the scalar type used by variables
template<Concepts::Variables V>
using ScalarType = LinearSystem::ScalarType<DofsType<V>>;

} // namespace Variables
} // namespace Dumux

#endif
