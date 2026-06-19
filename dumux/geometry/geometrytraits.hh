// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Geometry
 * \brief Traits telling the bounding box tree algorithms about a leaf geometry
 */
#ifndef DUMUX_GEOMETRY_GEOMETRY_TRAITS_HH
#define DUMUX_GEOMETRY_GEOMETRY_TRAITS_HH

#include <dune/geometry/axisalignedcubegeometry.hh>

namespace Dumux {

/*!
 * \ingroup Geometry
 * \brief Traits telling the bounding box tree algorithms about a leaf geometry
 *
 * The member \c axisAligned indicates whether a geometry coincides with its
 * axis-aligned bounding box. If so, the cheap bounding box test is already
 * exact and the more expensive primitive test (point-in-geometry,
 * geometry intersection) can be skipped during queries. Specialize this for
 * any geometry that is axis-aligned.
 */
template<class Geometry>
struct GeometryTraits
{ static constexpr bool axisAligned = false; };

//! \brief Specialization for Dune's axis-aligned cube geometry
template<class ct, int dim, int dimworld>
struct GeometryTraits<Dune::AxisAlignedCubeGeometry<ct, dim, dimworld>>
{ static constexpr bool axisAligned = true; };

} // end namespace Dumux

#endif
