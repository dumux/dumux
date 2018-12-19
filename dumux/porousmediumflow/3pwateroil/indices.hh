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
 * \ingroup ThreePWaterOilModel
 * \brief Defines the indices required for the 3p2cni model.
 */

#ifndef DUMUX_3P2CNI_INDICES_HH
#define DUMUX_3P2CNI_INDICES_HH

namespace Dumux {

/*!
 * \ingroup ThreePWaterOilModel
 * \brief The indices for the isothermal 3p2cni model.
 */
class ThreePWaterOilIndices
{
public:
    // present phases (-> 'pseudo' primary variable)
    static const int threePhases = 1; //!< All three phases are present
    static const int wPhaseOnly = 2; //!< Only the water phase is present
    static const int gnPhaseOnly = 3; //!< Only gas and NAPL phases are present
    static const int wnPhaseOnly = 4; //!< Only water and NAPL phases are present
    static const int gPhaseOnly = 5; //!< Only gas phase is present
    static const int wgPhaseOnly = 6; //!< Only water and gas phases are present

    // Primary variable indices
    static const int pressureIdx = 0; //!< Index for phase pressure in a solution vector
    static const int switch1Idx = 1; //!< Index 1 of saturation, mole fraction or temperature
    static const int switch2Idx = 2; //!< Index 2 of saturation, mole fraction or temperature

    // equation indices
    static const int conti0EqIdx = 0; //!< Index of the mass conservation equation for the water component
};

} // end namespace Dumux

#endif
