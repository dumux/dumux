// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FaceCenteredStaggeredDiscretization
 * \copydoc Dumux::ConsistentlyOrientedGrid
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_CONSISTENTLY_ORIENTED_GRID_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_CONSISTENTLY_ORIENTED_GRID_HH

#include <type_traits>

// forward declare
namespace Dune {
template<int dim, class Coordinates>
class YaspGrid;

template <int dim, class HostGrid, bool mapIndexStorage>
class SubGrid;

}

namespace Dumux {

/*!
 * \brief Helper type to determine whether a grid is guaranteed to be oriented consistently.
 *        This means that the intersection indices always correspond to the ones of a reference element
 *        or, in other words, the elements are never rotated.
 */
template<class T>
struct ConsistentlyOrientedGrid : public std::false_type {};

template<int dim, class Coords>
struct ConsistentlyOrientedGrid<Dune::YaspGrid<dim, Coords>> : public std::true_type {};

template<int dim, class Coords, bool mapIndexStorage>
struct ConsistentlyOrientedGrid<Dune::SubGrid<dim, Dune::YaspGrid<dim, Coords>, mapIndexStorage>> : public std::true_type {};


} // end namespace Dumux

#endif
