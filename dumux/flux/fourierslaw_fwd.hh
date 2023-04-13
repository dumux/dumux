// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Flux
 * \brief Fourier's law specialized for different discretization schemes
 *        This file contains the data which is required to calculate
 *        diffusive mass fluxes due to molecular diffusion with Fourier's law.
 */
#ifndef DUMUX_FLUX_FOURIERS_LAW_FWD_HH
#define DUMUX_FLUX_FOURIERS_LAW_FWD_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

// declaration of primary template
template <class TypeTag, class DiscretizationMethod>
class FouriersLawImplementation;

/*!
 * \ingroup Flux
 * \brief Evaluates the heat conduction flux according to Fouriers's law
 */
template <class TypeTag>
using FouriersLaw = FouriersLawImplementation<TypeTag, typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod>;

} // end namespace Dumux

#endif
