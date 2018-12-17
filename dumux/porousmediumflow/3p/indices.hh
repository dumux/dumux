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
 * \ingroup ThreePModel
 * \brief Defines the indices for the three-phase model.
 */

#ifndef DUMUX_3P_INDICES_HH
#define DUMUX_3P_INDICES_HH

namespace Dumux {

/*!
 * \ingroup ThreePModel
 * \brief The common indices for the isothermal three-phase model.
 */
class ThreePIndices
{
public:
    // Primary variable indices
    static constexpr int pressureIdx = 0; //!< index for gas phase pressure in a solution vector
    static constexpr int swIdx = 1; //!< index of water (more wetting than the other liquid) saturation
    static constexpr int snIdx = 2; //!< index of (e.g.) NAPL saturation

    // equation indices
    static constexpr int conti0EqIdx = 0; //!< index of first balance equation
};

} // end namespace Dumux

#endif
