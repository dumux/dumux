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
 * \brief Quantities required by the one-phase fully implicit model defined on a vertex.
 */

#ifndef DUMUX_1P_INITIAL_FLUIDSTATE_HH
#define DUMUX_1P_INITIAL_FLUIDSTATE_HH

namespace Dumux{
template<class Scalar, class FluidSystem, class FluidState, class EnergyVolVars>
class InitialFluidState
{
public:
    /*!
     * \brief complete the fluid state from an initial state.
     *
     * \param initFS
     * \return FluidState
     */
    static FluidState completeFluidState(const FluidState& initFS)
    {
        FluidState fluidState(initFS);

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(fluidState, /*phaseIdx=*/0);

        Scalar value = FluidSystem::density(fluidState, paramCache, /*phaseIdx=*/0);
        fluidState.setDensity(/*phaseIdx=*/0, value);

        value = FluidSystem::viscosity(fluidState, paramCache, /*phaseIdx=*/0);
        fluidState.setViscosity(/*phaseIdx=*/0, value);

        // compute and set the enthalpy
        value = EnergyVolVars::enthalpy(fluidState, paramCache, /*phaseIdx=*/0);
        fluidState.setEnthalpy(/*phaseIdx=*/0, value);

        return fluidState;
    }
};
}

#endif
