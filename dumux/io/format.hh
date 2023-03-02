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

#if __has_include(<format>) // cppcheck-suppress preprocessorErrorDirective
#include <format>
#endif

#include <dumux/io/format/fmt/format.h>
#include <dumux/io/format/fmt/ranges.h>

//! Formatting tools in the style of std::format (C++20)
namespace Dumux::Fmt {

#if __cpp_lib_format
// use std::format from C++20
using std::format;
using std::format_to;
using std::format_to_n;
using std::formatted_size;
using std::vformat;
using std::vformat_to;
using std::make_format_args;
#else
// use fallback fmt library
using Dumux::Detail::fmt::format;
using Dumux::Detail::fmt::format_to;
using Dumux::Detail::fmt::format_to_n;
using Dumux::Detail::fmt::formatted_size;
using Dumux::Detail::fmt::vformat;
using Dumux::Detail::fmt::vformat_to;
using Dumux::Detail::fmt::make_format_args;
#endif

} // end namespace Dumux::Fmt

#endif
