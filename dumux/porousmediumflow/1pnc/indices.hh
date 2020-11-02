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
 * \ingroup OnePNCModel
 * \brief Defines the primary variable and equation indices used by
 *        the 1pnc model
 */

#ifndef DUMUX_1PNC_INDICES_HH
#define DUMUX_1PNC_INDICES_HH

namespace Dumux {

/*!
 * \ingroup OnePNCModel
 * \brief The indices for the isothermal one-phase n-component model.
 *
 * \tparam phaseIdx The index of the fluid phase in the fluid system
 */
struct OnePNCIndices
{
    //! Reference index for mass conservation equation.
    static constexpr int conti0EqIdx = 0;
    //! Index for wetting/nonwetting phase pressure (depending on formulation) in a solution vector
    static constexpr int pressureIdx = 0;
};

} // end namespace Dumux

#endif
