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
 * \ingroup InputOutput
 * \brief Formatting based on the fmt-library which implements std::format of C++20
 *
 * For documentation of the functions, see https://en.cppreference.com/w/cpp/utility/format
 * For a documentation of formatting styles,
 * see https://en.cppreference.com/w/cpp/utility/format/formatter#Standard_format_specification
 *
 * Once std::format / C++20 is available, we can use the standard library here.
 */
#ifndef DUMUX_IO_FORMAT_HH
#define DUMUX_IO_FORMAT_HH

#include <dumux/io/format/fmt/format.h>
#include <dumux/io/format/fmt/ranges.h>

//! Formatting tools in the style of std::format (C++20)
namespace Dumux::Fmt {

using fmt::format;
using fmt::format_to;
using fmt::format_to_n;
using fmt::formatted_size;

} // end namespace Dumux::Fmt

#endif
