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
 * \ingroup Adaptive
 * \brief A free function for h-adaptivity.
 */
#ifndef DUMUX_ADAPTIVE_ADAPT_HH
#define DUMUX_ADAPTIVE_ADAPT_HH

#include "griddatatransfer.hh"

namespace Dumux {

/*!
 * \ingroup Adaptive
 * \brief Adapt the grid and reconstruct the user data
 *
 * \param grid The grid to adapt
 * \param dataTransfer A class that performs the data
 *                     transfer from the old to the new grid.
 */
template<class Grid>
bool adapt(Grid& grid, GridDataTransfer& dataTransfer)
{
    // Do pre-adaption step of the grid
    const bool mightCoarsen = grid.preAdapt();

    // Let the helper do storage of variables
    dataTransfer.store();

    // adapt the grid
    const bool refine = grid.adapt();

    // (Re-)construct variables to new grid
    dataTransfer.reconstruct();

    // delete markers in grid
    grid.postAdapt();

    // return boolean if grid has been changed
    return mightCoarsen || refine;
}

} // end namespace Dumux

#endif
