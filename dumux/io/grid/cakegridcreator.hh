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
 * \brief Provides a grid creator for a piece of cake grid
 */
#ifndef DUMUX_CAKE_GRID_CREATOR_HH
#define DUMUX_CAKE_GRID_CREATOR_HH

#warning "This header is deprecated. Use the new cakegridmanager."

#include "cakegridmanager.hh"

namespace Dumux {

/*!
 * \ingroup InputOutput
 * \brief Provides a grid creator with a method for creating creating vectors
 *        with polar Coordinates and one for creating a cartesian grid from
 *        these polar coordinates.
 */
template<class Grid>
using CakeGridCreator [[deprecated("Use CakeGridManager instead!")]] = CakeGridManager<Grid>;
} // end namespace Dumux

#endif
