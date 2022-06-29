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
 * \ingroup ConstraintSolvers
 * \brief Contains update solvers for secondary variables
 */

#ifndef DUMUX_CONSTRAINT_SOLVERS_MOBILITY_HH
#define DUMUX_CONSTRAINT_SOLVERS_MOBILITY_HH

namespace Dumux {

template < class MobilityType, class FluidState, class FluidMatrixInteraction>
void updateMobility(MobilityType mobility,
                    const FluidState& fluidState,
                    const FluidMatrixInteraction& fluidMatrixInteraction)
{
    const int wPhaseIdx = fluidState.wettingPhase();
    const int nPhaseIdx = 1 - wPhaseIdx;

    mobility[wPhaseIdx] = fluidMatrixInteraction.krw(fluidState.saturation(wPhaseIdx))
                        / fluidState.viscosity(wPhaseIdx);

    mobility[nPhaseIdx] = fluidMatrixInteraction.krn(fluidState.saturation(wPhaseIdx))
                        / fluidState.viscosity(nPhaseIdx);
}

template <class MobilityType, class FluidState, class FluidMatrixInteraction>
void updateMobilityMP(MobilityType mobility,
                      const FluidState& fluidState,
                      const FluidMatrixInteraction& fluidMatrixInteraction,
                      const int& wPhaseIdx,
                      const int& nPhaseIdx,
                      const int& numFluidPhases)
{
    const auto sw = fluidState.saturation(wPhaseIdx);
    const auto sn = fluidState.saturation(nPhaseIdx);

    // update the mobilities
    for (int phaseIdx = 0; phaseIdx < numFluidPhases; ++phaseIdx)
    {
        mobility[phaseIdx] = fluidMatrixInteraction.kr(phaseIdx, sw, sn)
                            / fluidState.viscosity(phaseIdx);
    }

}

template <class MobilityType, class RelativePermeabilityType, class FluidState, class FluidMatrixInteraction>
void updateMobilityMPNC(MobilityType mobility,
                        RelativePermeabilityType relativePermeability,
                        const FluidState& fluidState,
                        const FluidMatrixInteraction& fluidMatrixInteraction,
                        const int& wPhaseIdx,
                        const int& numFluidPhases)
{
    // update the relative permeabilities
    const auto relPerm = fluidMatrixInteraction.relativePermeabilities(fluidState, wPhaseIdx);
    std::copy(relPerm.begin(), relPerm.end(), relativePermeability.begin());

    // update the mobilities
    for (int phaseIdx = 0; phaseIdx < numFluidPhases; ++phaseIdx)
        mobility[phaseIdx] = relativePermeability[phaseIdx] / fluidState.viscosity(phaseIdx);
}

} // end namespace Dumux

#endif
