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
 * \brief Types and helpers for defining dynamic/static container sizes.
 */
#ifndef DUMUX_COMMON_SIZE_HH
#define DUMUX_COMMON_SIZE_HH

#include <concepts>
#include <type_traits>

namespace Dumux {
namespace Size {

/*!
 * \ingroup Common
 * \brief Represents dynamic container sizes.
 *        Can be used in contexts were both reserved and dynamic container sizes
 *        are supported. Any arithmetic operation to "compute" a container size
 *        from multiple values yields a dynamic size again.
 */
struct Dynamic
{
    friend constexpr auto operator<=>(const Dynamic&, const Dynamic&) = default;
    friend constexpr Dynamic operator+(const Dynamic&, const Dynamic&) { return {}; }
    friend constexpr Dynamic operator*(const Dynamic&, const Dynamic&) { return {}; }
    friend constexpr Dynamic operator+(const Dynamic&, const std::integral auto&) { return {}; }
    friend constexpr Dynamic operator*(const Dynamic&, const std::integral auto&) { return {}; }
    friend constexpr Dynamic operator+(const std::integral auto&, const Dynamic&) { return {}; }
    friend constexpr Dynamic operator*(const std::integral auto&, const Dynamic&) { return {}; }
};

inline constexpr Dynamic dynamic{};

} // namespace Size

namespace Concepts {

template<typename T>
concept Size = std::integral<T> or std::is_same_v<T, Dumux::Size::Dynamic>;

} // namespace Dumux::Concepts

namespace Size {

/*!
 * \ingroup Common
 * \brief Returns true if the given size `a` is able to represent the size `b`.
 */
inline constexpr bool isRepresentableBy(const Dynamic& a, const Dynamic& b) { return true; }
inline constexpr bool isRepresentableBy(const Dynamic& a, const std::integral auto& b) { return true; }
inline constexpr bool isRepresentableBy(const std::integral auto& a, const Dynamic& b) { return false; }
inline constexpr bool isRepresentableBy(const std::integral auto& a, const std::integral auto& b) { return a >= b; }

template<Concepts::Size auto v>
inline constexpr bool isDynamic = std::is_same_v<std::decay_t<decltype(v)>, Size::Dynamic>;

#ifndef DOXYGEN
namespace Detail {

template<Concepts::Size auto size,
         Concepts::Size auto exponent,
         Concepts::Size auto size0>
constexpr auto pow()
{
    if constexpr (isDynamic<size> || isDynamic<exponent> || isDynamic<size0>)
        return Size::dynamic;
    else
    {
        static_assert(exponent >= 0);
        if constexpr (exponent == 1)
            return size;
        else
            return pow<size*size0, exponent-1, size0>();
    }
}

} // namespace Detail
#endif // DOXYGEN

template<Concepts::Size auto size,
         Concepts::Size auto exponent>
constexpr auto pow()
{ return Detail::pow<size, exponent, size>(); }

} // namespace Size
} // namespace Dumux

#endif
