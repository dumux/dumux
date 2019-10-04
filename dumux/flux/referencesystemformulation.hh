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
 * \brief The reference frameworks and formulations available for splitting total fluxes into a advective and diffusive part
 */

#ifndef DUMUX_FLUX_REFERENCESYSTEMFORMULATION_HH
#define DUMUX_FLUX_REFERENCESYSTEMFORMULATION_HH

namespace Dumux {

/*!
 * \brief The formulations available for Fick's law related to the reference system
 * \ingroup Flux
 * \note The total flux of a component can be split into an advective and a diffusive part.
 *       In our framework, the advective part is based on a momentum balance. Standard momentum balances, e.g.,
 *       the Navier-Stokes equations or Darcy's law yield mass-averaged velocities (see Multicomponent Mass Transfer, Taylor & Krishna, 1993 \cite taylor1993a),
 *       therefore we use the appropriate formulation of Fick's law (mass averaged formulation) per default.
 *
 * This means that the diffusive fluxes are calculated with the mass fraction gradients and the unit of the fluxes is kg/s.
 * It is also possible to use a molar-averaged reference system, which can be beneficial, e.g.,
 * when it is known that the molar-averaged advective velocity would be zero. When using a molar-averaged reference velocity,
 * Fick's law is formulated with mole fraction gradients and the unit of the flux is moles/s. This means that depending on the reference system,
 * the units of the fluxes need to be adapted to be used in mass or mole balances.
*/
enum class ReferenceSystemFormulation
{
    massAveraged, molarAveraged
};

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
typename VolumeVariables::PrimaryVariables::value_type
massOrMoleFraction(const VolumeVariables& volVars, ReferenceSystemFormulation referenceSys, const int phaseIdx, const int compIdx)
{
    return (referenceSys == ReferenceSystemFormulation::massAveraged) ? volVars.massFraction(phaseIdx, compIdx) : volVars.moleFraction(phaseIdx, compIdx);
}

} // end namespace Dumux

#endif
