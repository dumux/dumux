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
 * \brief Basic traits that provide customization points for user types.
 */
#ifndef DUMUX_COMMON_TRAITS_HH
#define DUMUX_COMMON_TRAITS_HH

#include <type_traits>
#include <concepts>
#include <limits>

namespace Dumux::Traits {

//! Trait to allow registering custom types as arithmetics
template<typename T>
struct IsArithmetic : public std::false_type {};

//! Trait to allow registering custom types as views
template<typename T>
struct IsView : public std::false_type {};

//! Trait to define a value of type `T` that can be used to denote "undefined" values
template<typename T>
struct UndefinedValue;

template<std::floating_point T>
struct UndefinedValue<T>
: public std::integral_constant<T, std::numeric_limits<T>::signaling_NaN()>
{};

template<std::integral T>
struct UndefinedValue<T>
: public std::integral_constant<T, std::numeric_limits<T>::max()>
{};

} // namespace Dumux::Traits

#endif
