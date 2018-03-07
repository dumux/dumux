// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup ThreePModel
 * \brief Defines the indices for the three-phase model.
 */
#ifndef DUMUX_3P_INDICES_HH
#define DUMUX_3P_INDICES_HH

#include <dumux/common/properties.hh>

namespace Dumux
{

/*!
 * \ingroup ThreePModel
 * \brief The common indices for the isothermal three-phase model.
 *
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <class FluidSystem, int PVOffset = 0>
class ThreePIndices
{
public:
    // Phase indices
    static const int wPhaseIdx = FluidSystem::wPhaseIdx; //!< index of the wetting liquid phase
    static const int nPhaseIdx = FluidSystem::nPhaseIdx; //!< index of the nonwetting liquid phase
    static const int gPhaseIdx = FluidSystem::gPhaseIdx; //!< index of the gas phase


    // Primary variable indices
    static const int pressureIdx = PVOffset + 0; //!< index for gas phase pressure in a solution vector
    static const int swIdx = PVOffset + 1; //!< index of water (more wetting than the other liquid) saturation
    static const int snIdx = PVOffset + 2; //!< index of (e.g.) NAPL saturation


    // equation indices
    static const int conti0EqIdx = PVOffset    + wPhaseIdx; //!< index of the mass conservation equation for the water component
    static const int conti1EqIdx = conti0EqIdx + nPhaseIdx; //!< index of the mass conservation equation for the contaminant component
    static const int conti2EqIdx = conti0EqIdx + gPhaseIdx; //!< index of the mass conservation equation for the air component

    static const int contiWEqIdx = conti0EqIdx + wPhaseIdx; //!< index of the mass conservation equation for the water component
    static const int contiNEqIdx = conti0EqIdx + nPhaseIdx; //!< index of the mass conservation equation for the contaminant component
    static const int contiGEqIdx = conti0EqIdx + gPhaseIdx; //!< index of the mass conservation equation for the air component
};

}

#endif
