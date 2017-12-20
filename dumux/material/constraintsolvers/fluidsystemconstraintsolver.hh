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
 * \ingroup ConstraintSolver
 * \brief Computes the composition of all phases of a N-phase,
 *        N-component fluid system assuming that all N phases are
 *        present. The composition is actually retrieved from a
 *        function in the fluidsystem.
 */
#ifndef DUMUX_MISCIBILITY_CONSTRAINT_SOLVER_HH
#define DUMUX_MISCIBILITY_CONSTRAINT_SOLVER_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/common/valgrind.hh>

namespace Dumux {
/*!
 * \ingroup ConstraintSolver
 * \brief Computes the composition of all phases from a function in the fluidsystem.
 *
 *        This is basically an interface in order to use ConstraintSolver with fluidsystems,
 *        which specify the solubilities / miscibilities / mole fractions in stead of
 *        some measure for chemical potential.
 *
 * The constraint solver assumes the following quantities to be set:
 *
 * - temperatures of *all* phases
 * - saturations of *all* phases
 * - pressures of *all* phases
 *
 *  After calling the solve() method the following quantities
 *  are calculated in addition:
 *
 * - density, molar density, molar volume of *all* phases
 * - composition in mole and mass fractions and molarities of *all* phases
 * - mean molar masses of *all* phases
 * - if the setViscosity parameter is true, also dynamic viscosities of *all* phases
 * - if the setEnthalpy parameter is true, also specific enthalpies of *all* phases
 */
template <class Scalar, class FluidSystem>
class FluidSystemConstraintSolver
{
    static constexpr int numPhases = FluidSystem::numPhases;

public:
    /*!
     * \brief @copybrief Dumux::FluidSystemConstraintSolver
     *
     * \param fluidState A container with the current (physical) state of the fluid
     * \param paramCache A container for iterative calculation of fluid composition
     * \param setViscosity Should the viscosity be set in the fluidstate?
     * \param setEnthalpy Should the enthalpy be set in the fluidstate?
     */
    template <class FluidState, class ParameterCache>
    static void solve(FluidState & fluidState,
                      ParameterCache & paramCache,
                      const bool setViscosity,
                      const bool setEnthalpy)
    {

        // In this function the actual work is done.
        // Either tables for solubility or functional relations are used therein.
        // One way or the other, this function needs to set all the mole fractions
        // in the fluidstate.
        FluidSystem::calculateEquilibriumMoleFractions(fluidState, paramCache);

        for (int phaseIdx=0; phaseIdx<numPhases; phaseIdx++){
            Scalar value = FluidSystem::density(fluidState, paramCache, phaseIdx);
            fluidState.setDensity(phaseIdx, value);

            if (setViscosity) {
                value = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
                fluidState.setViscosity(phaseIdx, value);
            }

            if (setEnthalpy) {
                value = FluidSystem::enthalpy(fluidState, paramCache, phaseIdx);
                fluidState.setEnthalpy(phaseIdx, value);
            }
        }
    }
};

} // end namespace Dumux

#endif
