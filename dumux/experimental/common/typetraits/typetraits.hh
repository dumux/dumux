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
#ifndef DUMUX_EXPERIMENTAL_TYPE_TRAITS_HH
#define DUMUX_EXPERIMENTAL_TYPE_TRAITS_HH

#include <type_traits>

namespace Dumux {

/*!
 * \brief Function that performs no operation.
 */
inline constexpr auto noop = [] (auto...) {};
using Noop = decltype(noop);


} // end namespace Dumux
#endif
