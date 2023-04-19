// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Flux
 * \brief Forchheimer's law specialized for different discretization schemes
 *        This file contains the data which is required to calculate
 *        volume and mass fluxes of fluid phases over a face of a finite volume by means
 *        of the Forchheimer approximation. Specializations are provided for the different discretization methods.
 */
#ifndef DUMUX_FLUX_FORCHHEIMERS_LAW_FWD_HH
#define DUMUX_FLUX_FORCHHEIMERS_LAW_FWD_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/flux/forchheimervelocity.hh>

namespace Dumux {

// definition of primary template
template <class TypeTag, class VelocityLaw, class DiscretizationMethod>
class ForchheimersLawImplementation
{
    static_assert(
        GetPropType<TypeTag, Properties::GridGeometry>::discMethod == DiscretizationMethods::cctpfa || 
        GetPropType<TypeTag, Properties::GridGeometry>::discMethod == DiscretizationMethods::box, 
        "Forchheimer only implemented for cctpfa or box!"
    );
};

/*!
 * \ingroup Flux
 * \brief Evaluates the normal component of the Forchheimer velocity on a (sub)control volume face.
 * \note Specializations are provided for the different discretization methods.
 * These specializations are found in the headers included below.
 */
template <class TypeTag>
using ForchheimersLaw = ForchheimersLawImplementation<
    TypeTag,
    ForchheimerVelocity<
        GetPropType<TypeTag, Properties::Scalar>,
        GetPropType<TypeTag, Properties::GridGeometry>,
        GetPropType<TypeTag, Properties::FluxVariables>
    >,
    typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod
>;

} // end namespace Dumux

#endif
