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
 * \brief Types and functions related to multi-dimensional indexing.
 */
#ifndef DUMUX_COMMON_MULTI_INDEX_HH
#define DUMUX_COMMON_MULTI_INDEX_HH

#include <tuple>
#include <array>
#include <ostream>
#include <concepts>
#include <type_traits>

namespace Dumux {
namespace Concepts {

//! Concept for indices defined at compile-time
template<typename T>
concept StaticIndex = requires {
    { T::get() };
    std::integral<std::decay_t<decltype(T::get())>>;
};


//! Concept for multi-indices
template<typename T>
concept MultiIndex = requires(const T& t) {
    { T::size() } -> std::integral;
    { t.template get<0>() };
    std::integral<std::decay_t<decltype(t.template get<0>())>>;
};


//! Concept for supported index types
template<typename T>
concept Index = std::integral<T> or StaticIndex<T> or MultiIndex<T>;

} // namespace Concepts


//! Default implementation of a compile-time index
template<std::size_t i>
struct StaticIndex
{
    static constexpr std::size_t get()
    { return i; }
};

namespace Indices {

//! Convenience variables for compile-time indices
template<std::size_t i>
inline constexpr StaticIndex<i> _i{};

inline constexpr auto _0 = _i<0>;
inline constexpr auto _1 = _i<1>;
inline constexpr auto _2 = _i<2>;
inline constexpr auto _3 = _i<3>;
inline constexpr auto _4 = _i<4>;
inline constexpr auto _5 = _i<5>;

} // namespace Indices


#ifndef DOXYGEN
namespace Detail {

template<std::integral I>
void printIndex(std::ostream& s, const I& index)
{ s << index; }

template<Concepts::StaticIndex I>
void printIndex(std::ostream& s, const I&)
{ s << I::get(); }

template<std::size_t i = 0, Concepts::MultiIndex I>
void printIndex(std::ostream& s, const I& index)
{
    if constexpr (i == 0)
        s << "(";

    printIndex(s, index.template get<i>());
    if constexpr (i < I::size() - 1)
    {
        s << ", ";
        printIndex<i+1>(s, index);
    }

    if constexpr (i == 0)
        s << ")";
}

} // namespace Detail
#endif // DOXYGEN


/*!
 * \ingroup Common
 * \brief Class that represents an index in a multi-dimensional data structure.
 */
template<Concepts::Index... Indices> requires(sizeof...(Indices) > 0)
class MultiIndex
{
    static constexpr std::size_t numIndices = sizeof...(Indices);

public:
    using IndexTuple = std::tuple<Indices...>;

    template<int i>
    using IndexType = std::tuple_element_t<i, IndexTuple>;

    constexpr MultiIndex(const Indices&... indices)
    : indices_(indices...)
    {}

    constexpr MultiIndex(const IndexTuple& indices)
    : indices_(indices)
    {}

    static constexpr std::size_t size()
    { return numIndices; }

    template<int i>
    constexpr decltype(auto) get() const
    { return std::get<i>(indices_); }

    template<std::integral I, int i>
    constexpr I getAs() const
    {
        if constexpr (Concepts::StaticIndex<IndexType<i>>)
            return static_cast<I>(get<i>().get());
        else if constexpr (std::integral<IndexType<i>>)
            return static_cast<I>(get<i>());
        else
            static_assert("Casting only possible for static or integral indices");
    }

    friend std::ostream& operator<<(std::ostream& s, const MultiIndex& i)
    {
        Detail::printIndex(s, i);
        return s;
    }

private:
    IndexTuple indices_;
};

template<Concepts::Index... Indices>
MultiIndex(const std::tuple<Indices...>& tuple) -> MultiIndex<Indices...>;


#ifndef DOXYGEN
namespace Detail {

template<int current, int start, int stop, Concepts::MultiIndex I, Concepts::Index... Indices>
inline constexpr auto getSubMultiIndex(const I& index, const std::tuple<Indices...>& indices)
{
    if constexpr (current == stop)
        return MultiIndex{
            std::tuple_cat(indices, std::make_tuple(index.template get<current>()))
        };
    else
        return getSubMultiIndex<current+1, start, stop>(
            index,
            std::tuple_cat(indices, std::make_tuple(index.template get<current>()))
        );
}

template<int current, int start, int stop, Concepts::MultiIndex I>
inline constexpr auto getSubMultiIndex(const I& index)
{
    if constexpr (current < start)
        return getSubMultiIndex<current+1, start, stop>(index);
    else
        return getSubMultiIndex<current+1, start, stop>(
            index,
            std::make_tuple(index.template get<current>())
        );
}

} // namespace Detail
#endif // DOXYGEN


//! Returns the sub-multi-index between the given index bounds
template<int start, int stop = -1, Concepts::MultiIndex I>
inline constexpr auto getSubMultiIndex(const I& index)
{
    static_assert(I::size() > 0);

    constexpr int actualStop = stop == -1 ? I::size() - 1 : stop;
    static_assert(actualStop < I::size());
    static_assert(start >= 0 && start <= actualStop);
    return Detail::getSubMultiIndex<0, start, stop>(index);
}


#ifndef DOXYGEN
namespace Detail {

template<int current, Concepts::MultiIndex MI, std::integral I>
void copyMultiIndex(const MI& index, std::array<I, MI::size()>& arr)
{
    static_assert(current < MI::size());
    arr[current] = index.template getAs<I, current>();
    if constexpr (current < MI::size() - 1)
        copyMultiIndex<current+1>(index, arr);
}

} // namespace Detail
#endif // DOXYGEN


//! Convert a multi index into an array
template<Concepts::MultiIndex I>
std::array<std::size_t, I::size()> asArray(const I& index)
{
    std::array<std::size_t, I::size()> result;
    Detail::copyMultiIndex<0>(index, result);
    return result;
}

} // namespace Dumux

#endif
