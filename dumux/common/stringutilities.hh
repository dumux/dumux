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
 * \brief Helpers for working with strings
 */

#ifndef DUMUX_COMMON_STRING_UTILITIES_HH
#define DUMUX_COMMON_STRING_UTILITIES_HH

#include <vector>
#include <string_view>
#include <tuple>

namespace Dumux {

/*
 * \brief Tokenize a string splitting at a given delimiter
 * \ingroup Common
 * \param str a string to be tokenized
 * \param delim a set of delimiter characters separating the tokens
 * \note this works without copying the original string, make sure that string
 *       does not go out of scope!
 * \note The delimiter characters cannot appear within any of the tokens
 *
 * Examples:
 *  - tokenize("bla&foo&&&bar", "&") -> {"bla", "foo", "bar"}
 *  - tokenize("    a  b   ", " ") -> {"a", "b"}
 *  - tokenize(str, " \n\t") -> split at whitespace (and removes all whitespace)
 */
std::vector<std::string_view> tokenize(std::string_view str, std::string_view delim)
{
    const auto token = [&](std::size_t pos){
        const auto start = str.find_first_not_of(delim, pos);
        const auto end = str.find_first_of(delim, start);
        return std::make_pair(start, end);
    };

    std::vector<std::string_view> tokens;
    for (auto [start, end] = token(0); start != end; std::tie(start, end) = token(end))
        tokens.emplace_back(str.substr(start, end-start));
    return tokens;
}

/*
 * \brief Split a string at a given seperator string
 * \ingroup Common
 * \param str a string to be split
 * \param delim a set of delimiter characters separating the tokens
 * \param removeEmpty if empty string between two separator should be removed from the result
 * \note this works without copying the original string, make sure that string
 *       does not go out of scope!
 *
 * Examples:
 *  - split("| Hello | world |", "|") -> {"", " Hello ", " world ", ""}
 *  - split("| Hello | world |", "|", true) -> {" Hello ", " world "}
 *  - split("blafoobla", "foo") -> {"bla", "bla"}
 *  - split("blafoobla", "bla") -> {"", "foo", ""}
 *  - split("blafoobla", "bla", true) -> {"foo"}
 */
std::vector<std::string_view> split(std::string_view str, std::string_view delim, bool removeEmpty = false)
{
    const auto token = [&](std::size_t pos){
        const auto end = str.find(delim, pos);
        return std::make_pair(pos, end);
    };

    std::vector<std::string_view> split;
    for (auto [start, end] = token(0); start != std::string::npos; std::tie(start, end) = token(start))
    {
        if (auto&& sub = str.substr(start, end-start); !removeEmpty || !sub.empty())
            split.emplace_back(std::move(sub));

        // skip to after the delimiter
        start = (end != std::string::npos) ? end + delim.size() : end;
    }
    return split;
}

} // end namespace Dumux

#endif
