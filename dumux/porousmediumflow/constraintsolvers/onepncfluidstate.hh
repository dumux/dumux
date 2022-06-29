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

#ifndef DUMUX_CONSTRAINT_SOLVERS_ONEPNCFLUIDSTATE_HH
#define DUMUX_CONSTRAINT_SOLVERS_ONEPNCFLUIDSTATE_HH

namespace Dumux {
/*!
*
* \brief Sets complete fluid state.
* \param elemSol A vector containing all primary variables connected to the element
* \param problem The object specifying the problem which ought to
*                be simulated
* \param element An element which contains part of the control volume
* \param scv The sub-control volume
* \param fluidState A container with the current (physical) state of the fluid
* \param solidState A container with the current (physical) state of the solid
*/
template<class Volvars, class ElemSol, class Problem, class Element, class Scv, class FluidState, class SolidState>
void completeFluidState(Volvars& volVars,
                        const ElemSol &elemSol,
                        const Problem& problem,
                        const Element& element,
                        const Scv &scv,
                        FluidState& fluidState,
                        SolidState& solidState,
                        const int& pressureIdx,
                        const int& numFluidComps,
                        const bool& useMoles = true)
{
    using Scalar = typename Volvars::Scalar;
    using FluidSystem = typename Volvars::FluidSystem;
    volVars.updateTemperature(elemSol, problem, element, scv, fluidState, solidState);
    fluidState.setSaturation(0, 1.0);

    const auto& priVars = elemSol[scv.localDofIndex()];
    fluidState.setPressure(0, priVars[pressureIdx]);

    // Set fluid state mole fractions
    if (useMoles)
    {
        Scalar sumMoleFracNotMainComp = 0;
        for (int compIdx = 1; compIdx < numFluidComps; ++compIdx)
        {
            fluidState.setMoleFraction(0, compIdx, priVars[compIdx]);
            sumMoleFracNotMainComp += priVars[compIdx];
        }
        fluidState.setMoleFraction(0, 0, 1.0 - sumMoleFracNotMainComp);
    }
    else
    {
        // for mass fractions we only have to set numComponents-1 mass fractions
        // the fluid state will internally compute the remaining mass fraction
        for (int compIdx = 1; compIdx < numFluidComps; ++compIdx)
            fluidState.setMassFraction(0, compIdx, priVars[compIdx]);
    }

    typename FluidSystem::ParameterCache paramCache;
    paramCache.updateAll(fluidState);

    Scalar rho = FluidSystem::density(fluidState, paramCache, 0);
    Scalar rhoMolar = FluidSystem::molarDensity(fluidState, paramCache, 0);
    Scalar mu = FluidSystem::viscosity(fluidState, paramCache, 0);

    fluidState.setDensity(0, rho);
    fluidState.setMolarDensity(0, rhoMolar);
    fluidState.setViscosity(0, mu);

    // compute and set the enthalpy
    Scalar h = volVars.enthalpy(fluidState, paramCache, 0);
    fluidState.setEnthalpy(0, h);
}
}
#endif
