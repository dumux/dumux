// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Flux
 * \brief Fick's law specilized for different discretization schemes.
 *        This file contains the data which is required to calculate
 *        diffusive mass fluxes due to molecular diffusion with Fick's law.
 */
#ifndef DUMUX_DISCRETIZATION_FICKS_LAW_HH
#define DUMUX_DISCRETIZATION_FICKS_LAW_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/flux/referencesystemformulation.hh>

namespace Dumux {

// forward declaration
template <class TypeTag, DiscretizationMethod discMethod, ReferenceSystemFormulation referenceSystem>
class FicksLawImplementation;

/*!
 * \ingroup Flux
 * \brief Evaluates the diffusive mass flux according to Fick's law
 */
template <class TypeTag, ReferenceSystemFormulation referenceSystem =  ReferenceSystemFormulation::massAveraged>
using FicksLaw = FicksLawImplementation<TypeTag, GetPropType<TypeTag, Properties::FVGridGeometry>::discMethod, referenceSystem>;

/*!
 * \ingroup Flux
 * \brief evaluates the density to be used in Fick's law based on the reference system
 */
template<class VolumeVariables>
typename VolumeVariables::PrimaryVariables::value_type
massOrMolarDensity(const VolumeVariables& volVars, ReferenceSystemFormulation referenceSys, const int phaseIdx)
{
    return (referenceSys == ReferenceSystemFormulation::massAveraged) ? volVars.density(phaseIdx) : volVars.molarDensity(phaseIdx);
}

/*!
 * \ingroup Flux
 * \brief returns the mass or mole fraction to be used in Fick's law based on the reference system
 */
template<class VolumeVariables>
typename VolumeVariables::PrimaryVariables::value_type massOrMoleFraction(const VolumeVariables& volVars, ReferenceSystemFormulation referenceSys, const int phaseIdx, const int compIdx)
{
    return (referenceSys == ReferenceSystemFormulation::massAveraged) ? volVars.massFraction(phaseIdx, compIdx) : volVars.moleFraction(phaseIdx, compIdx);
}

} // end namespace Dumux

#include <dumux/flux/cctpfa/fickslaw.hh>
#include <dumux/flux/ccmpfa/fickslaw.hh>
#include <dumux/flux/box/fickslaw.hh>
#include <dumux/flux/staggered/freeflow/fickslaw.hh>

#endif
