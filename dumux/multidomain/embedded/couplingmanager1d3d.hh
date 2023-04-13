// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedCoupling
 * \brief Coupling manager for low-dimensional domains embedded in the bulk domain.
 */

#ifndef DUMUX_MULTIDOMAIN_EMBEDDED_COUPLINGMANAGER_1D3D_HH
#define DUMUX_MULTIDOMAIN_EMBEDDED_COUPLINGMANAGER_1D3D_HH

namespace Dumux {

/*!
 * \ingroup EmbeddedCoupling
 * \brief Manages the coupling between bulk elements and lower dimensional elements
 *        Point sources on each integration point are computed by an AABB tree.
 */
template<class MDTraits, class CouplingMode>
class Embedded1d3dCouplingManager;

} // end namespace Dumux

#include <dumux/multidomain/embedded/couplingmanager1d3d_line.hh>
#include <dumux/multidomain/embedded/couplingmanager1d3d_average.hh>
#include <dumux/multidomain/embedded/couplingmanager1d3d_surface.hh>
#include <dumux/multidomain/embedded/couplingmanager1d3d_kernel.hh>
#include <dumux/multidomain/embedded/couplingmanager1d3d_projection.hh>

#endif
