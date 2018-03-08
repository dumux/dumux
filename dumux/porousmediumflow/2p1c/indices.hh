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
 * \ingroup TwoPOneCModel
 * \copydoc Dumux::TwoPOneCIndices
 *
 */
#ifndef DUMUX_2P1C_INDICES_HH
#define DUMUX_2P1C_INDICES_HH

#include <dumux/common/properties.hh>

namespace Dumux {

/*!
 * \ingroup TwoPOneCModel
 * \brief The indices for the two-phase one-component model.
 *
 * \tparam FluidSystem The fluid system class
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <class FluidSystem, int PVOffset = 0>
class TwoPOneCIndices
{
public:
    // Phase indices
    static const int wPhaseIdx = FluidSystem::wPhaseIdx; //!< Index of the wetting phase.
    static const int nPhaseIdx = FluidSystem::nPhaseIdx; //!< Index of the non-wetting phase.

    // Present phases (-> 'pseudo' primary variable)
    static const int twoPhases = 1; //!< Both wetting and non-wetting phase are present.
    static const int wPhaseOnly = 2; //!< Only the wetting phase is present.
    static const int nPhaseOnly = 3; //!< Only non-wetting phase is present.

    // Primary variable indices
    static const int pressureIdx = PVOffset + 0; //!< Index for phase pressure in a solution vector.
    static const int switch1Idx = PVOffset + 1; //!< Index of saturation or temperature.

    // Equation indices
    static const int conti0EqIdx = PVOffset    + 0; //!< Index of the mass conservation equation for the water component.
    static const int energyEqIdx = PVOffset + 1; //<! The index for energy in equation vectors.
};

}

#endif
