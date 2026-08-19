// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Core
 * \brief Helpers for deprecation
 */

#ifndef DUMUX_COMMON_DEPRECATED_HH
#define DUMUX_COMMON_DEPRECATED_HH

#include <utility>

#include <dune/common/ftraits.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/std/type_traits.hh>
#include <dumux/common/numeqvector.hh>

namespace Dumux {

#ifndef DOXYGEN // hide from doxygen
// Helper classes/functions for deprecation
// Each implementation has to state after which release
// it will be removed. Implementations in the Deprecated
// namespace will be removed without
// deprecation after their usage in the code expired,
// so most likely you don't want to use this in your code
namespace Deprecated {

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif // __clang__

//! Get the grid discretization from a problem or an element discretization
//! (e.g. FVElementGeometry, FEElementDiscretization), preferring the new
//! gridDiscretization() interface where available and otherwise falling
//! back to the gridGeometry() interface.
template<class T>
decltype(auto) gridGeometry(const T& t)
{
    if constexpr (requires { t.gridDiscretization(); })
        return t.gridDiscretization();
    else
        return t.gridGeometry();
}

#ifdef __clang__
#pragma clang diagnostic pop
#endif  // __clang__


} // end namespace Deprecated
#endif

} // end namespace Dune
#endif
