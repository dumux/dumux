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
 * \ingroup MPNCModel
 * \brief The primary variable and equation indices for the MpNc model.
 */

#ifndef DUMUX_MPNC_INDICES_HH
#define DUMUX_MPNC_INDICES_HH

namespace Dumux {

/*!
 * \ingroup MPNCModel
 * \brief The primary variable and equation indices for the MpNc model.
 *
 * \tparam FluidSystem The fluid system class
 * \tparam numEqBalance Number of balance equations: all transport equations and the constraint equations
 */
template <int numPhases, int numEqBalance>
struct MPNCIndices
{
    /*!
     * \brief Index of the saturation of the first phase in a vector
     *        of primary variables.
     *
     * \note The following (numPhases - 1) primary variables represent the
     *       saturations for the phases [1, ..., numPhases - 1]
     */
    static const unsigned int s0Idx = numEqBalance - numPhases;

    /*!
     * \brief Index of the first phase' pressure in a vector of primary variables.
     */
    static const unsigned int p0Idx = numEqBalance  - 1;

    /*!
     * \brief Index of the first phase NCP equation.
     * \note The index for the remaining phases are consecutive.
     */
    static const unsigned int phase0NcpIdx =  numEqBalance - numPhases;

    static const unsigned int fug0Idx = 0;
    static const unsigned int conti0EqIdx = 0;
    static const unsigned int moleFrac00Idx = 0;
};

}

#endif
