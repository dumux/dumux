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
 * \brief  Defines the indices for the one-phase fully implicit model.
 */
#ifndef DUMUX_1P_MIMETIC_INDICES_HH
#define DUMUX_1P_MIMETIC_INDICES_HH

namespace Dumux
{
// \{
/*!
 * \ingroup NavierStokesModel
 * \ingroup ImplicitIndices
 * \brief The common indices for the isothermal stokes model.
 *
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <class TypeTag, int PVOffset = 0>
struct OnePMimeticIndices
{
    static const int conti0EqIdx = PVOffset + 0; //index for the mass balance
    static const int pressureIdx = PVOffset + 0; //index of the primary variable

    static constexpr int faceOffset = GET_PROP_VALUE(TypeTag, NumEqCellCenter); //!< Index of the momentum balance equation

    static constexpr int faceFluxBalanceIdx = PVOffset + faceOffset; //!< Index of the momentum balance equation
    static constexpr int facePressureIdx = faceFluxBalanceIdx;

    static constexpr int phaseIdx = 0; //!< Index of the fluid phase (required to use the same fluid system in coupled models)
};

// \}
} // end namespace

#endif
