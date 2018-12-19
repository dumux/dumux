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
 * \ingroup RichardsModel
 * \brief Index names for the Richards model.
 */

#ifndef DUMUX_RICHARDS_INDICES_HH
#define DUMUX_RICHARDS_INDICES_HH

namespace Dumux {

/*!
 * \ingroup RichardsModel
 * \brief Index names for the Richards model.
 */

struct RichardsIndices
{
    //! Primary variable index for the wetting phase pressure
    static constexpr int pressureIdx = 0;
    static constexpr int switchIdx = 0;

    //! Equation index for the mass conservation of the wetting phase
    static constexpr int conti0EqIdx = 0;

    // present phases (-> 'pseudo' primary variable)
    static constexpr int liquidPhaseOnly = 1; //!< Only the liquid phase is present
    static constexpr int gasPhaseOnly = 2; //!< Only the gas phase is present
    static constexpr int bothPhases = 3; //!< Both phases are present
};

} // end namespace Dumux

#endif
