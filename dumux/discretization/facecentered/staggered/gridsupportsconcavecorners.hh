// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FaceCenteredStaggeredDiscretization
 * \copydoc Dumux::GridSupportsConcaveCorners
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_GRID_SUPPORTS_CONCAVE_CORNERS_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_GRID_SUPPORTS_CONCAVE_CORNERS_HH

#include <type_traits>

// forward declare
namespace Dune {

template<int dim, class Coordinates>
class YaspGrid;

}

namespace Dumux {

/*!
 * \brief Type trait to determine if a grid supports concave corners (e.g. by cutting out a hole from the domain interior)
 */
template<class T>
struct GridSupportsConcaveCorners : public std::true_type {};

template<int dim, class Coords>
struct GridSupportsConcaveCorners<Dune::YaspGrid<dim, Coords>> : public std::false_type {};

} // end namespace Dumux

#endif
