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
 * \brief  A grid manager for dune-subgrid.
 */
#ifndef DUMUX_SUBGRID_MANAGER_HH
#define DUMUX_SUBGRID_MANAGER_HH

#if HAVE_DUNE_SUBGRID
#warning "This header is deprecated and will be removed after release 3.1. Use gridmanager_sub.hh"

#include <dumux/io/grid/gridmanager_sub.hh>

namespace Dumux {


/*!
 * \ingroup InputOutput
 * \brief A grid manager for dune-subgrid
 */
template<class HostGrid, class HostGridManager = GridManager<HostGrid>>
using SubgridManager [[deprecated("Use GridManager<SubGrid> instead. Will be removed after release 3.1")]]
 = GridManager<Dune::SubGrid<HostGrid::dimension, HostGrid>>;

} // end namespace Dumux

#endif // HAVE_DUNE_SUBGRID
#endif // DUMUX_SUBGRID_MANAGER_HH
