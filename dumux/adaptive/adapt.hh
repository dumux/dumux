// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
bool adapt(Grid& grid, GridDataTransfer<Grid>& dataTransfer)
{
    // Do pre-adaption step of the grid
    const bool mightCoarsen = grid.preAdapt();

    // Let the helper do storage of variables
    dataTransfer.store(grid);

    // adapt the grid
    const bool refine = grid.adapt();

    // (Re-)construct variables to new grid
    dataTransfer.reconstruct(grid);

    // delete markers in grid
    grid.postAdapt();

    // return boolean if grid has been changed
    return mightCoarsen || refine;
}

} // end namespace Dumux

#endif
