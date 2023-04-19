// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup InputOutput
 * \brief Vtk output precision options available in Dumux
 */
#ifndef DUMUX_IO_VTK_PRECISION_HH
#define DUMUX_IO_VTK_PRECISION_HH

#include <array>
#include <string_view>
#include <dune/grid/io/file/vtk/common.hh>

namespace Dumux::Vtk {

using Dune::VTK::Precision;

/*!
 * \ingroup InputOutput
 * \brief Maps a string (e.g. from input) to a Dune precision type
 */
inline Precision stringToPrecision(std::string_view precisionName)
{
    // this should really be constexpr but GCC <= 7.2 has a bug which
    // doesn't allow string_view to be constexpr
    static const std::array<std::pair<std::string_view, Precision>, 5> nameToPrecision
    {{
        { "Float32", Precision::float32 },
        { "Float64", Precision::float64 },
        { "UInt32", Precision::uint32 },
        { "UInt8", Precision::uint8 },
        { "Int32", Precision::int32 },
    }};

    for (const auto& [name, precision] : nameToPrecision)
        if (name == precisionName)
            return precision;

    DUNE_THROW(Dune::InvalidStateException, "Unknown precision type " << precisionName);
}

} // end namespace Dumux::Vtk

#endif
