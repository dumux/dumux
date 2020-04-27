// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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

#ifndef DUMUX_ICP_CO2TABLES_LABORATORY_HH
#define DUMUX_ICP_CO2TABLES_LABORATORY_HH

// ## The CO2 tables (`co2tableslaboratory.hh`)
//
// This file contains the __co2table class__ which forwards to tabulated properties of CO2 according to Span and Wagner 1996.
// The real work (creating the tables) is done by some external program by Span and Wagner 1996 which provides the ready-to-use tables.
//
// [[content]]
//
// [[codeblock]]

#include <assert.h>
namespace Dumux::ICP {
#include "co2valueslaboratory.inc"
}// end namespace Dumux::ICP
// [[/codeblock]]
// [[/content]]
#endif
