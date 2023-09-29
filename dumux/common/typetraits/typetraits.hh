// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Typetraits
 * \brief Type traits.
 */
#ifndef DUMUX_TYPE_TRAITS_HH
#define DUMUX_TYPE_TRAITS_HH

#include <type_traits>

namespace Dumux {

/*!
 * \brief Template which always yields a false value
 * \tparam T Some type.
 */
template<typename T>
struct AlwaysFalse : public std::false_type {};

/*!
 * \brief Function that performs no operation.
 */
inline constexpr auto noop = [] (auto...) {};
using Noop = decltype(noop);

/*!
 * \brief Helper template to select type T if it is not void
 *        or fall back to the given default type otherwise.
 */
template<typename Default, typename T>
using NonVoidOr = std::conditional_t<!std::is_void_v<T>, T, Default>;

} // end namespace Dumux

#endif
