/*****************************************************************************
 *   Copyright (C) 2009-2010 by Melanie Darcis                               *
 *   Copyright (C) 2010 by Andreas Lauser                                    *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/**
 * \file
 *
 * \brief Provides the class with the tabulated values of CO2 for the
 *        benchmark3 problem
 */
#ifndef DUMUX_HETEROGENEOUS_CO2TABLES_HH
#define DUMUX_HETEROGENEOUS_CO2TABLES_HH

#include <cassert>
#include <dumux/material/components/co2tablereader.hh>

namespace Dumux
{
namespace HeterogeneousCO2Tables
{
// the real work is done by some external program which provides
// ready-to-use tables.
#include "co2values.inc"
}
}

#endif
