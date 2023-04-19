// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Flux
 * \brief Fick's law specialized for different discretization schemes.
 *        This file contains the data which is required to calculate
 *        diffusive mass fluxes due to molecular diffusion with Fick's law.
 */
#ifndef DUMUX_FLUX_FICKS_LAW_FWD_HH
#define DUMUX_FLUX_FICKS_LAW_FWD_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>

#include <dumux/flux/referencesystemformulation.hh>

namespace Dumux {

// declaration of primary template
template <class TypeTag, class DiscretizationMethod, ReferenceSystemFormulation referenceSystem>
class FicksLawImplementation;

/*!
 * \ingroup Flux
 * \brief Evaluates the diffusive mass flux according to Fick's law
 */
template <class TypeTag, ReferenceSystemFormulation referenceSystem =  ReferenceSystemFormulation::massAveraged>
using FicksLaw = FicksLawImplementation<TypeTag, typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod, referenceSystem>;

} // end namespace Dumux

#endif
