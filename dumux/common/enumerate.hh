// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Core
 * \brief A Python-like enumerate function
 */

#ifndef DUMUX_COMMON_ENUMERATE_HH
#define DUMUX_COMMON_ENUMERATE_HH

#include <tuple>

namespace Dumux {

/*!
 * \brief A Python-like enumerate function
 * \param iterable Some iterable type with begin/end (e.g. std::vector)
 * \note From Nathan Reed (CC BY 4.0): http://reedbeta.com/blog/python-like-enumerate-in-cpp17/
 *
 * Usage example: for (const auto& [i, item] : enumerate(list))
 */
template <typename Range,
          typename RangeIterator = decltype(std::begin(std::declval<Range>())),
          typename = decltype(std::end(std::declval<Range>()))>
constexpr auto enumerate(Range&& iterable)
{
    struct Iterator
    {
        std::size_t i;
        RangeIterator iter;
        bool operator != (const Iterator & other) const { return iter != other.iter; }
        void operator ++ () { ++i; ++iter; }
        auto operator * () const { return std::tie(i, *iter); }
    };

    struct Iterable
    {
        Range iterable;
        auto begin() { return Iterator{ 0, std::begin(iterable) }; }
        auto end() { return Iterator{ 0, std::end(iterable) }; }
    };

    return Iterable{ std::forward<Range>(iterable) };
}

} // end namespace Dumux

#endif
