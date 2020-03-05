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
 * \ingroup TypeTraits
 * \brief Utilities for template meta programming
 */
#ifndef DUMUX_COMMON_TYPETRAITS_UTILITY_HH
#define DUMUX_COMMON_TYPETRAITS_UTILITY_HH

#include <cstddef>
#include <utility>

namespace Dumux {

/*
 * \ingroup TypeTraits
 * \brief create a variadic template from indexed types
 * \tparam V a variadic template that we want to create
 * \tparam T an indexed type (type that gets an index as template parameter)
 * \tparam U the list of indices
 */
template <template<typename... Args> class Variadic, template<std::size_t> class Indexed, class U>
struct makeFromIndexedType;

template <template<typename... Args> class Variadic, template<std::size_t> class Indexed, std::size_t... IndexSeq>
struct makeFromIndexedType<Variadic, Indexed, std::index_sequence<IndexSeq...>>
{
    using type = Variadic<Indexed<IndexSeq>...>;
};

namespace Detail {
    template <class Seq1, std::size_t offset, class Seq2> struct ConcatSeq;

    template <std::size_t ... Is1, std::size_t offset, std::size_t ... Is2>
    struct ConcatSeq<std::index_sequence<Is1...>, offset, std::index_sequence<Is2...>>
    {
        using type = std::index_sequence<Is1..., (offset + Is2)...>;
    };
}

/*
 * \ingroup TypeTraits
 * \brief create an integer sequence from 0 to n-1, omitting one specific number e
 * \tparam n number of integers in complete sequence before omitting
 * \tparam e value of integer to be omitted
 *
 * example: makeIncompleteIntegerSequence<3, 1> = [0, 2]
 * example: makeIncompleteIntegerSequence<4, 4> = [0, 1, 2, 3]
 *
 * see https://stackoverflow.com/questions/27124920/compile-time-generate-integer-sequence-with-one-left-out for details
 */
template <std::size_t n, std::size_t e>
using makeIncompleteIntegerSequence =
        typename Detail::ConcatSeq<decltype(std::make_index_sequence<e>{}), e + 1, decltype(std::make_index_sequence<(n > e) ? (n - e - 1) : 0>{})>::type;

/*
 * \ingroup TypeTraits
 * \brief add an offset to an index sequence
 * \tparam offset the offset
 * \tparam is the index sequence
 *
 * see https://stackoverflow.com/a/35625414
 */
template <std::size_t offset, std::size_t ... is>
constexpr std::index_sequence<(offset + is)...> addOffsetToIndexSequence(std::index_sequence<is...>)
{ return {}; }

/*
 * \ingroup TypeTraits
 * \brief create an index sequence starting from an offset
 * \tparam offset the offset
 * \tparam n the length of the sequence
 *
 * example: makeIndexSequenceWithOffset<2, 3> = [2,3,4]
 *
 * see https://stackoverflow.com/a/35625414
 */
template <std::size_t offset, std::size_t n>
constexpr auto makeIndexSequenceWithOffset()
{
    return addOffsetToIndexSequence<offset>(std::make_index_sequence<n>{});
}

} // end namespace Dumux

#endif
