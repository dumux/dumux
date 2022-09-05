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
 * \brief Strategies for mapping of multi indices.
 */
#ifndef DUMUX_COMMON_INDEX_STRATEGIES_HH
#define DUMUX_COMMON_INDEX_STRATEGIES_HH

#include <array>
#include <numeric>
#include <functional>
#include <type_traits>

#include <dumux/experimental/new_assembly/dumux/common/multiindex.hh>

namespace Dumux {
namespace Concepts {

//! Concept for the block sizes of multi-dimensional arrays.
template<typename T>
concept BlockSizes = MultiIndex<T>;


//! Concpt for index strategies, converting a `MultiIndex` via `operator[]`
template<typename T, typename I>
concept IndexStrategy = MultiIndex<I> and requires(const T& t, const I& i) {
    { t[i] };
    Concepts::MultiIndex<std::decay_t<decltype(t[i])>>;
};

} // namespace Concepts


//! A type usable to define block sizes.
template<Concepts::Index... Sizes>
using BlockSizes = MultiIndex<Sizes...>;


//! Index strategy for blocked multi-dimensional arrays
template<int blockDepth = -1>
class BlockedIndexStrategy
{
    static constexpr bool hasArbitraryBlockDepth = blockDepth == -1;

public:
    template<Concepts::MultiIndex I> requires(
        hasArbitraryBlockDepth or I::size() == blockDepth)
    constexpr const I& operator[](const I& i) const
    { return i; }
};


#ifndef DOXYGEN
namespace Detail {

template<int blockDepth, std::integral I, std::size_t size>
std::size_t getBlockOffset(const std::array<I, size>& sizes)
{
    static_assert(blockDepth < size);
    return std::accumulate(
        sizes.begin() + blockDepth,
        sizes.end(),
        1,
        std::multiplies{}
    );
}

template<int blockDepth, Concepts::MultiIndex MI, std::integral I>
std::size_t getFlatIndex(const MI& index, const std::array<I, MI::size()>& blockSizes)
{
    static_assert(blockDepth < MI::size());
    if constexpr (blockDepth == MI::size() - 1)
        return index.template getAs<std::size_t, blockDepth>();
    else
        return index.template getAs<std::size_t, blockDepth>()
               * getBlockOffset<blockDepth+1>(blockSizes)
               + getFlatIndex<blockDepth+1>(index, blockSizes);
}

} // namespace Detail
#endif // DOXYGEN


//! Index strategy for flat access to multi-dimensional data
template<int blockDepth>
class FlatIndexStrategy
{
public:
    template<Concepts::Index... Sizes> requires(
        sizeof...(Sizes) == blockDepth)
    FlatIndexStrategy(const Sizes&... sizes)
    : sizes_(asArray(MultiIndex{sizes...}))
    {}

    template<Concepts::BlockSizes Sizes> requires(
        Sizes::size() == blockDepth)
    FlatIndexStrategy(const Sizes& sizes)
    : sizes_(asArray(sizes))
    {}

    template<Concepts::MultiIndex I> requires(
        I::size() == blockDepth)
    constexpr auto operator[](const I& i) const
    { return MultiIndex{Detail::getFlatIndex<0>(i, sizes_)}; }

private:
    std::array<std::size_t, blockDepth> sizes_;
};

template<Concepts::BlockSizes Sizes>
FlatIndexStrategy(const Sizes& sizes) -> FlatIndexStrategy<Sizes::size()>;

} // namespace Dumux

#endif
