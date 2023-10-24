// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Core
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
 * \ingroup Core
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
// cppcheck-suppress internalAstError
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
