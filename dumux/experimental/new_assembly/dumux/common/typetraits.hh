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
 * \brief Type traits.
 */
#ifndef DUMUX_COMMON_TYPE_TRAITS_HH
#define DUMUX_COMMON_TYPE_TRAITS_HH

#include <type_traits>

#include <dumux/experimental/new_assembly/dumux/common/traits.hh>

namespace Dumux {

#ifndef DOXYGEN
namespace Detail {

template<typename T, std::size_t s = sizeof(T)>
std::false_type isIncomplete(T*);
std::true_type isIncomplete(...);

} // namespace Detail
#endif // DOXYGEN


template<typename T>
inline constexpr bool isIncomplete = decltype(Detail::isIncomplete(std::declval<T*>()))::value;


template<typename T>
inline constexpr bool isComplete = !isIncomplete<T>;


template<typename T, typename I = std::size_t>
inline constexpr bool isIndexable = requires(T& t, I& i) { { t[i] }; };


//! The value type obtained from indexing
template<typename T> requires isIndexable<T>
using IndexedType = std::decay_t<decltype(std::declval<T>()[0])>;


#ifndef DOXYGEN
namespace Detail {

template<typename T>
struct IndexDepth
: public std::integral_constant<std::size_t, 0>
{};

template<typename T> requires isIndexable<T>
struct IndexDepth<T>
: public std::integral_constant<
    std::size_t,
    1 + IndexDepth<IndexedType<T>>::value
> {};

} // namespace Detail
#endif // DOXYGEN

template<typename T>
inline constexpr std::size_t indexDepth = Detail::IndexDepth<T>::value;


//! A value of type `T` that can be used to denote undefined values
template<typename T> requires(isComplete<Traits::UndefinedValue<T>>)
inline constexpr T undefined = Traits::UndefinedValue<T>::value;


//! Type that can be used to indicate that a type should be selected automatically
struct Auto {};

//! Type that can be used to indicate that something should by dynamically sized
struct DynamicSize
{
    friend constexpr auto operator<=>(const DynamicSize&, const DynamicSize&) = default;
    friend constexpr DynamicSize operator+(const DynamicSize&, const DynamicSize&) { return {}; }
    friend constexpr DynamicSize operator*(const DynamicSize&, const DynamicSize&) { return {}; }
    friend constexpr DynamicSize operator+(const DynamicSize&, const std::integral auto&) { return {}; }
    friend constexpr DynamicSize operator*(const DynamicSize&, const std::integral auto&) { return {}; }
    friend constexpr DynamicSize operator+(const std::integral auto&, const DynamicSize&) { return {}; }
    friend constexpr DynamicSize operator*(const std::integral auto&, const DynamicSize&) { return {}; }
};

namespace Size {

template<auto v>
inline constexpr bool isDynamic = std::is_same_v<std::decay_t<decltype(v)>, DynamicSize>;
inline constexpr DynamicSize dynamic{};

template<auto v, auto n, auto v0 = v>
constexpr auto pow()
{
    using Size::isDynamic;
    if constexpr (isDynamic<v> || isDynamic<n> || isDynamic<v0>)
        return DynamicSize{};
    else
    {
        static_assert(n >= 0);
        if constexpr(n <= 1)
            return v;
        else
            return pow<v*v0, n-1, v>();
    }
}

} // namespace Size
} // namespace Dumux

#endif
