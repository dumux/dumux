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
