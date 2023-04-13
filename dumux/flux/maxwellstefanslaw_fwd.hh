// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Flux
 * \brief This file contains the data which is required to calculate
 *        diffusive mass fluxes due to molecular diffusion with Maxwell-Stefan's law.
 */
#ifndef DUMUX_FLUX_MAXWELL_STEFAN_LAW_FWD_HH
#define DUMUX_FLUX_MAXWELL_STEFAN_LAW_FWD_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/flux/referencesystemformulation.hh>

namespace Dumux {

// forward declaration
template <class TypeTag, class DiscretizationMethod, ReferenceSystemFormulation referenceSystem>
class MaxwellStefansLawImplementation;

/*!
 * \ingroup Flux
 * \brief Evaluates the diffusive mass flux according to Maxwell Stefan's law
 */
template <class TypeTag, ReferenceSystemFormulation referenceSystem =  ReferenceSystemFormulation::massAveraged>
using MaxwellStefansLaw = MaxwellStefansLawImplementation<TypeTag, typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod, referenceSystem>;

} // end namespace Dumux

#endif
