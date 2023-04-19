// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
