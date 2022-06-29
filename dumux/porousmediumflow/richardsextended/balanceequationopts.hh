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
 * \brief Traits class to set options used by the local residual when
 *        when evaluating the balance equations.
 */
#ifndef DUMUX_POROUSMEDIUMFLOW_RICHARDS_BALANCE_EQUATION_OPTIONS_HH
#define DUMUX_POROUSMEDIUMFLOW_RICHARDS_BALANCE_EQUATION_OPTIONS_HH

namespace Dumux {

/*!
 * \ingroup RichardsModel
 * \brief Traits class to set options used by the local residual when
 *        when evaluating the balance equations.
 */
template <class FluidSystem>
struct RichardsBalanceEquationOptions
{
    /*
     * The main component in the liquid phase (phase 0) is always balanced
     * and the main component in the gas phase (phase 1) is never balanced
     * For the index convention see RichardsVolumeVariables.
     */
    static constexpr bool mainComponentIsBalanced(int phaseIdx)
    { return phaseIdx == FluidSystem::phase0Idx; }
};

} // end namespace Dumux

#endif
