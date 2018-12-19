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
 * \ingroup ThreePThreeCModel
 * \brief Defines the indices required for the three-phase three-component fully implicit model.
 */

#ifndef DUMUX_3P3C_INDICES_HH
#define DUMUX_3P3C_INDICES_HH

namespace Dumux {

/*!
 * \ingroup ThreePThreeCModel
 * \brief The indices for the isothermal three-phase three-component model.
 */
class ThreePThreeCIndices
{
public:
    // present phases (-> 'pseudo' primary variable)
    static constexpr int threePhases = 1; //!< All three phases are present
    static constexpr int wPhaseOnly = 2; //!< Only the water phase is present
    static constexpr int gnPhaseOnly = 3; //!< Only gas and NAPL phases are present
    static constexpr int wnPhaseOnly = 4; //!< Only water and NAPL phases are present
    static constexpr int gPhaseOnly = 5; //!< Only gas phase is present
    static constexpr int wgPhaseOnly = 6; //!< Only water and gas phases are present

    // Primary variable indices
    static constexpr int pressureIdx = 0; //!< index for gas phase pressure in a solution vector
    static constexpr int switch1Idx = 1; //!< index 1 of saturation or mole fraction
    static constexpr int switch2Idx = 2; //!< index 2 of saturation or mole fraction

    //! index for gas phase pressure in a solution vector
    static constexpr int pgIdx = pressureIdx;
    //! index of either the saturation of the wetting phase or the mole fraction secondary component if a phase is not present
    static constexpr int sOrX1Idx = switch1Idx;
    //! index of either the saturation of the nonwetting phase or the mole fraction secondary component if a phase is not present
    static constexpr int sOrX2Idx = switch2Idx;

    // equation indices
    static constexpr int conti0EqIdx = 0;
};

} // end namespace Dumux

#endif
