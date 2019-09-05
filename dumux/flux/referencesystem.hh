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

#ifndef DUMUX_FLUX_REFERENCESYSTEM_HH
#define DUMUX_FLUX_REFERENCESYSTEM_HH

namespace Dumux {

    /*!
     * \brief The formulations available for Fick's law related to the reference system
     * \ingroup Flux
     * \note The total flux of a component can be split into an advective and a diffusive part. The advective part is in our framework determined by a momentum balance. Standard momentum balances, e.g. Navier-Stokes and Darcy's law give back mass-averaged velocites (see Multicomponent Mass Transfer, R. Taylor u. R. Krishna. J.), therefore as default we use the formulation of Fick's law which matches to these velocities (mass averaged formulation).
     *
     * This means that the diffusive fluxes are calculated with the mass fraction gradients and the unit of the fluxes is in kg/s. It is also possible to use a molar averaged reference system, which can be benefitial e.g. when it is known that the molar averaged advective velocity would be zero. When using a molar averaged reference velocity Fick's law is formulated with mole fraction gradients and the unit of the flux is mol/s. This means that depending on the reference system the units of the fluxes need to be adapted to be used in mass or mole balances.
     */
    enum class ReferenceSystemFormulation
    {
        massAveraged, molarAveraged
    };

} // end namespace Dumux

#endif
