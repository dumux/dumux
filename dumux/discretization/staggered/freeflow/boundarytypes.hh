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
 * \ingroup StaggeredDiscretization
 * \copydoc Dumux::StaggeredFreeFlowBoundaryTypes
 */
#ifndef STAGGERED_FREEFLOW_BOUNDARY_TYPES_HH
#define STAGGERED_FREEFLOW_BOUNDARY_TYPES_HH

#warning "This header is deprecated. Use dumux/freeflow/navierstokes/boundarytypes.hh"

#include <dumux/freeflow/navierstokes/boundarytypes.hh>

namespace Dumux {

template <int numEq>
using StaggeredFreeFlowBoundaryTypes [[deprecated("Use NavierStokesBoundaryTypes instead. Will be removed after 3.3")]] = NavierStokesBoundaryTypes<numEq>;

} // end namespace Dumux

#endif
