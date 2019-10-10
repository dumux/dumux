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
 * \ingroup Discretization
 * \brief Free function to get the local view of a grid cache object
 */

#ifndef DUMUX_LOCAL_VIEW_HH
#define DUMUX_LOCAL_VIEW_HH

namespace Dumux {

/*!
 * \ingroup Discretization
 * \brief Free function to get the local view of a grid cache object
 * \note A local object is only functional after calling its bind/bindElement method.
 * \tparam GridCache the grid caching type (such as GridGeometry)
 * \param gridCache the grid caching object we want to localView from
 */
template<class GridCache>
inline typename GridCache::LocalView localView(const GridCache& gridCache)
{ return typename GridCache::LocalView(gridCache); }

} // end namespace Dumux

#endif
