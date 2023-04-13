// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Geometry
 * \brief Compute the center point of a convex polytope geometry or a random-access container of corner points
 */
#ifndef DUMUX_GEOMETRY_CENTER_HH
#define DUMUX_GEOMETRY_CENTER_HH

#include <numeric>

namespace Dumux {

/*!
 * \ingroup Geometry
 * \brief The center of a given list of corners
 */
template<class Corners>
typename Corners::value_type center(const Corners& corners)
{
    using Pos = typename Corners::value_type;
    auto center = std::accumulate(corners.begin(), corners.end(), Pos(0.0));
    center /= corners.size();
    return center;
}

} // end namespace Dumux

#endif
