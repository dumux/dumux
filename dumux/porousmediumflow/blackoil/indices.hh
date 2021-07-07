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
 * \ingroup BlackOilModel
 * \brief Defines the indices required for the three-phase three-component fully implicit model.
 */

#ifndef DUMUX_BLACK_OIL_INDICES_HH
#define DUMUX_BLACK_OIL_INDICES_HH

namespace Dumux {

/*!
 * \ingroup BlackOilModel
 * \brief The primary variable and equation indices for the black-oil model.
 */
class BlackOilIndices
{
public:


    // Primary variable indices
    static constexpr int pressureIdx = 0; //!< index for gas phase pressure in a solution vector
    static constexpr int switch1Idx = 1; //!< index 1 of saturation
    static constexpr int switch2Idx = 2; //!< index 2 of saturation

    // Alternative names for the indices
    static constexpr int saturation1Idx = switch1Idx;
    static constexpr int saturation2Idx = switch2Idx;


    // equation indices
    static constexpr int conti0EqIdx = 0;
};

} // end namespace Dumux

#endif
