// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Flux
 * \brief Hooke's law specialized for different discretization schemes.
 *        This computes the stress tensor and surface forces resulting from mechanical deformation.
 */
#ifndef DUMUX_DISCRETIZATION_HOOKES_LAW_FWD_HH
#define DUMUX_DISCRETIZATION_HOOKES_LAW_FWD_HH

#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup Flux
 * \brief This computes the stress tensor and surface forces resulting from mechanical deformation.
 * \note Specializations are provided for the different discretization methods.
 * These specializations are found in the headers included below.
 */
template <class Scalar, class GridGeometry, class DiscretizationMethod = typename GridGeometry::DiscretizationMethod>
class HookesLaw;

} // end namespace Dumux

#endif
