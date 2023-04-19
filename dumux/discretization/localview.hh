// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
