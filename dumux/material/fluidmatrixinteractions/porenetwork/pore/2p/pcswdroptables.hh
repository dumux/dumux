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
/**
  * \file
  * \brief Provides the class with the tabulated values of pc-sw curve for droplets
  */

#ifndef DUMUX_PCSWTABLES_HH
#define DUMUX_PCSWTABLES_HH

#include "pcswtablereader.hh"

namespace Dumux {

//Provides the class with the tabulated values of pc-sw curve
namespace PcSwTables {

#include "pcswtables.inc"
} //end namespace pcswtables
} //end namespace Dumux
#endif
