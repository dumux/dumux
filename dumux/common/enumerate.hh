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

    struct Iteratable
    {
        Range iterable;
        auto begin() { return Iterator{ 0, std::begin(iterable) }; }
        auto end() { return Iterator{ 0, std::end(iterable) }; }
    };

    return Iteratable{ std::forward<Range>(iterable) };
}

} // end namespace Dumux

#endif
