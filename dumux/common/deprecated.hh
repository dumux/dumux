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
#include <ranges>

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

#ifdef __clang__
#pragma clang diagnostic pop
#endif  // __clang__

// Helper concept detecting new interface.
// Remove after release 3.11
template<class GG, class GI>
concept hasRangeOfPeriodicallyMappedDofs = requires (GG gg, GI gi) {
    { gg.periodicallyMappedDofs(gi) } -> std::ranges::range;
};

// helper function triggering deprecation warning if grid geometry does not implement new interface or program uses old interface.
// Remove after release 3.11
template<class Dof>
[[deprecated("The periodicDofMap with values of single periodic dofs will not be supported after release 3.11. Use/define in custom grid geometry a periodicDofMap with a range of dofs as value type, as well as accessor periodicallyMappedDofs returing a std::views::all/single of the storage if possible")]]
inline const std::ranges::range auto wrapSinglePeriodicDof(const Dof& dof)
{
    return std::ranges::views::single(dof);
}

// Helper function to access values from map of periodic dofs.
// Remove after release 3.11
template<typename T>
inline const std::ranges::range auto ensureRangeOfPeriodicDofs(const T& t);

template<std::ranges::range T>
inline const std::ranges::range auto ensureRangeOfPeriodicDofs(const T& range)
{ return range; }

template<typename T> requires (!std::ranges::range<T>)
inline const std::ranges::range auto ensureRangeOfPeriodicDofs(const T& entry)
{ return wrapSinglePeriodicDof(entry); }

// Helper functions to access ranges of periodically mapped dofs through common interface.
// Remove after release 3.11
template<class GG, class GI>
inline const std::ranges::range auto rangeOfPeriodicallyMappedDofs(const GG& gg, const GI& gi)
{
    if constexpr (hasRangeOfPeriodicallyMappedDofs<GG, GI>)
        return ensureRangeOfPeriodicDofs(gg.periodicallyMappedDofs(gi));
    else
        return wrapSinglePeriodicDof(gg.periodicallyMappedDof(gi));
}

} // end namespace Deprecated
#endif

} // end namespace Dune
#endif
