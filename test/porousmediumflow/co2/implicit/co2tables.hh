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
/**
 * \file
 * \ingroup CO2Tests
 * \brief Provides the class with the tabulated values of CO2 density
 *        and enthalpy.
 */

#ifndef DUMUX_HETEROGENEOUS_CO2TABLES_HH
#define DUMUX_HETEROGENEOUS_CO2TABLES_HH

#include <dumux/material/components/co2tablereader.hh>

namespace Dumux {

// Provides the class with the tabulated values of CO2 density and enthalpy
namespace HeterogeneousCO2Tables {

#ifndef DOXYGEN // hide from doxygen
// the real work is done by some external program which provides
// ready-to-use tables.
#include "co2values.inc"
#endif

} // end namespace HeterogeneousCO2Tables
} // end namespace Dumux

#endif
