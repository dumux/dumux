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
 * \brief Common concepts
 */
#ifndef DUMUX_COMMON_CONCEPTS_HH
#define DUMUX_COMMON_CONCEPTS_HH

#include <concepts>
#include <ranges>

namespace Dumux::Concepts {

//! Concept for a range with a specified range value type
template<typename T, typename ValueType>
concept RangeOf = std::ranges::range<T> and std::same_as<std::ranges::range_value_t<T>, ValueType>;

} // namespace Dumux::Concepts

#endif
