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
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesIndices
 */
#ifndef DUMUX_NAVIERSTOKES_MASS_1P_INDICES_HH
#define DUMUX_NAVIERSTOKES_MASS_1P_INDICES_HH

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief The common indices for the isothermal Navier-Stokes mass conservation model.
 */
struct NavierStokesMassOnePIndices
{
    static constexpr int conti0EqIdx = 0; //!< Index of the first (total for pure-fluid systems) mass balance equation
    static constexpr int pressureIdx = conti0EqIdx; //!< Index of the pressure
    static constexpr int phasefield1EqIdx = 1;
    static constexpr int phasefield2EqIdx = 2;
    static constexpr int phasefield3EqIdx = 3;
    static constexpr int phi1Idx = phasefield1EqIdx;
    static constexpr int phi2Idx = phasefield2EqIdx;
    static constexpr int phi3Idx = phasefield3EqIdx;
    static constexpr int u1TransportEqIdx = 4;
    static constexpr int u2TransportEqIdx = 5;
    static constexpr int u3TransportEqIdx = 6;
    static constexpr int u1Idx = u1TransportEqIdx;
    static constexpr int u2Idx = u2TransportEqIdx;
    static constexpr int u3Idx = u3TransportEqIdx;
};

} // end namespace Dumux

#endif
