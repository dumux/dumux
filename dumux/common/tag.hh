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
 * \brief Helper class to create (named and comparable) tagged types
 */
#ifndef DUMUX_COMMON_TAG_HH
#define DUMUX_COMMON_TAG_HH

#include <sstream>
#include <ostream>
#include <type_traits>
#include <dune/common/classname.hh>
#include <dumux/common/typetraits/isvalid.hh>

namespace Dumux::Utility {

/*!
 * \ingroup Common
 * \brief Helper class to create (named and comparable) tagged types
 * Tags any given type. The tagged type is equality comparable and can be written to streams.
 * A custom name can be provided by implementing the `name()` member function.
 */
template<class T>
struct Tag {};

//! Tags are equality comparable and return true if the tagged types are equal
template<class T1, class T2>
inline constexpr bool operator==(Tag<T1>, Tag<T2>)
{ return std::is_same_v<T1, T2>; }

template<class T1, class T2>
inline constexpr bool operator!=(Tag<T1>, Tag<T2>)
{ return !std::is_same_v<T1, T2>; }

namespace Detail {
constexpr auto hasName = isValid([](auto&& t) -> decltype(t.name(), void()) {});
} // end namespace Detail

//! Return the class name of the tagged type calling t.name()
template<class T, std::enable_if_t<std::is_base_of_v<Tag<T>, T>, int> = 0>
auto operator<<(std::ostream& os, const T& t)
-> std::enable_if_t<decltype(Detail::hasName(t))::value, std::ostream&>
{ os << t.name(); return os; }

//! Return the class name of the tagged type calling Dune::className if t.name() doesn't exist
template<class T, std::enable_if_t<std::is_base_of_v<Tag<T>, T>, int> = 0>
auto operator<<(std::ostream& os, const T& t)
-> std::enable_if_t<!decltype(Detail::hasName(t))::value, std::ostream&>
{
    const auto fullName = Dune::className<T>();

    // strip all namespace qualifiers
    const auto pos = fullName.rfind("::");
    const auto name = pos != std::string::npos ? fullName.substr(pos+2) : fullName;

    os << name;
    return os;
}

} // end namespace Dumux::Utility

#endif
