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
 * \ingroup OnePModel
 * \brief  Defines the indices for the one-phase fully implicit model.
 */

#ifndef DUMUX_1P_INDICES_HH
#define DUMUX_1P_INDICES_HH

namespace Dumux {
// \{

/*!
 * \ingroup OnePModel
 * \brief Indices for the one-phase model.
 *
 * \tparam offset The first index in a primary variable vector.
 */
template<int offset = 0>
struct OnePIndices
{
    static const int PVOffset = offset;      //!< the first index in primary variable vector
    static const int conti0EqIdx = PVOffset; //!< index for the mass balance
    static const int pressureIdx = PVOffset; //!< index of the primary variable
};

// \}
} // end namespace

#endif
