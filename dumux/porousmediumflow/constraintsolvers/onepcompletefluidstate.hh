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

#ifndef DUMUX_CONSTRAINT_SOLVERS_ONEPFLUIDSTATE_HH
#define DUMUX_CONSTRAINT_SOLVERS_ONEPFLUIDSTATE_HH

namespace Dumux {
/*!
 * \brief Sets complete fluid state
 *
 * \param elemSol A vector containing all primary variables connected to the element
 * \param problem The object specifying the problem which ought to
 *                be simulated
 * \param element An element which contains part of the control volume
 * \param scv The sub-control volume
 * \param fluidState A container with the current (physical) state of the fluid
 * \param solidState A container with the current (physical) state of the solid
 */
template<class VolVars, class ElemSol, class Problem, class Element, class Scv, class FluidState, class SolidState>
void completeFluidState(VolVars& volVars,
                        const ElemSol& elemSol,
                        const Problem& problem,
                        const Element& element,
                        const Scv& scv,
                        FluidState& fluidState,
                        SolidState& solidState)
{
    volVars.updateTemperature(elemSol, problem, element, scv, fluidState, solidState);
    fluidState.setSaturation(/*phaseIdx=*/0, 1.);

    const auto& priVars = elemSol[scv.localDofIndex()];
    static constexpr int pressureIdx = VolVars::Indices::pressureIdx;
    fluidState.setPressure(/*phaseIdx=*/0, priVars[pressureIdx]);

    // saturation in a single phase is always 1 and thus redundant
    // to set. But since we use the fluid state shared by the
    // immiscible multi-phase models, so we have to set it here...
    fluidState.setSaturation(/*phaseIdx=*/0, 1.0);

    using FluidSystem = typename VolVars::FluidSystem;
    typename FluidSystem::ParameterCache paramCache;
    paramCache.updatePhase(fluidState, /*phaseIdx=*/0);

    double value = FluidSystem::density(fluidState, paramCache, /*phaseIdx=*/0);
    fluidState.setDensity(/*phaseIdx=*/0, value);

    value = FluidSystem::viscosity(fluidState, paramCache, /*phaseIdx=*/0);
    fluidState.setViscosity(/*phaseIdx=*/0, value);

    // compute and set the enthalpy
    value = volVars.enthalpy(fluidState, paramCache, /*phaseIdx=*/0);
    fluidState.setEnthalpy(/*phaseIdx=*/0, value);
    }
}
#endif
