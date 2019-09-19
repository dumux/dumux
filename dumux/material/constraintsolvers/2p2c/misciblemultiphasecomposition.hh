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
 * \brief Computes the composition of all phases of a two-phase,
 *        two-component fluid system assuming that all two phases are
 *        present and ideal mixtures (meaning fugacity coefficients
 *        are independent from composition)
 */
#ifndef DUMUX_2P2C_MISCIBLE_MULTIPHASE_COMPOSITION_HH
#define DUMUX_2P2C_MISCIBLE_MULTIPHASE_COMPOSITION_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/common/valgrind.hh>

namespace Dumux {
/*!
 * \ingroup ConstraintSolvers
 * \brief Computes the composition of all phases of a two-phase,
 *        two-component fluid system assuming that all two phases are
 *        present and ideal mixtures (meaning fugacity coefficients
 *        are independent from composition)
 *
 * The constraint solver assumes the following quantities to be set:
 *
 * - temperatures of *all* phases
 * - pressures of *all* phases
 *
 * It also assumes that the mole/mass fractions of all phases sum up
 * to 1 and that the fugacities of each component are equal in both phases.
 * After calling the solve() method the following quantities
 * are calculated in total:
 *
 * - composition in mole fractions of *all* phases
 *   \f$x^\kappa_\alpha\f$, \f$x^\kappa_\beta\f$
 * - fugacity coefficients of *all* components in *all* phases
 *   \f$\Phi^\kappa_\alpha\f$, \f$\Phi^\kappa_\beta\f$
 * - density and molar density of *all* phases
 *   \f$\rho_\alpha\f$, \f$\rho_\beta\f$, \f$\rho_{mol, \alpha}\f$, \f$\rho_{mol, \beta}\f$,
 *
 * Setting those quantities also guarantees by design of the fluid state
 * that the following quantities are set:
 *
 * - molar volume of *all* phases
 *   \f$V_{mol, \alpha}\f$
 * - mean molar masses of *all* phases \f$M_\alpha\f$
 * - mass fractions and molarities of *all* phases
 *   \f$X^\kappa_\alpha\f$, \f$c^\kappa_\alpha\f$
 */
template <class Scalar, class FluidSystem>
class TwoPTwoCMiscibleMultiPhaseComposition
{
    static constexpr int numPhases = FluidSystem::numPhases;
    static constexpr int numComponents = FluidSystem::numComponents;
    static const int numMajorComponents = FluidSystem::numPhases;

    enum
    {
        comp0Idx = FluidSystem::comp0Idx,
        comp1Idx = FluidSystem::comp1Idx,
        phase0Idx = FluidSystem::phase0Idx,
        phase1Idx = FluidSystem::phase1Idx
    };

public:
    /*!
     * \brief @copybrief Dumux::MiscibleMultiPhaseComposition
     *
     * \param fluidState A container with the current (physical) state of the fluid
     * \param paramCache A container for iterative calculation of fluid composition
     * \param knownPhaseIdx The index of the phase with known properties
     */
    template <class FluidState, class ParameterCache>
    static void solve(FluidState &fluidState,
                      ParameterCache &paramCache)
    {
#ifndef NDEBUG
        // this solver can only handle fluid systems which
        // assume ideal mixtures of all fluids
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            assert(FluidSystem::isIdealMixture(phaseIdx));

        }
#endif

        // compute all fugacity coefficients
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            paramCache.updatePhase(fluidState, phaseIdx);

            // since we assume ideal mixtures, the fugacity
            // coefficients of the components cannot depend on
            // composition, i.e. the parameters in the cache are valid
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                Scalar fugCoeff = FluidSystem::fugacityCoefficient(fluidState, paramCache, phaseIdx, compIdx);
                fluidState.setFugacityCoefficient(phaseIdx, compIdx, fugCoeff);
            }
        }

        // we define mole equilibrium ratios to make the following
        // calculations more readable
        Scalar k10 = fluidState.fugacityCoefficient(phase0Idx, comp0Idx)*fluidState.pressure(phase0Idx)/
                     (fluidState.fugacityCoefficient(phase1Idx, comp0Idx)*fluidState.pressure(phase1Idx));
        Scalar k11 = fluidState.fugacityCoefficient(phase0Idx, comp1Idx)*fluidState.pressure(phase0Idx)/
                     (fluidState.fugacityCoefficient(phase1Idx, comp1Idx)*fluidState.pressure(phase1Idx));

        // using the assumptions above, the equilibrium mole fractions are calculated
        fluidState.setMoleFraction(phase0Idx, comp0Idx, (1 - k11)/(k10 - k11));
        fluidState.setMoleFraction(phase1Idx, comp0Idx, fluidState.moleFraction(phase0Idx,comp0Idx) * k10);
        fluidState.setMoleFraction(phase0Idx, comp1Idx, 1.0 - fluidState.moleFraction(phase0Idx,comp0Idx));
        fluidState.setMoleFraction(phase1Idx, comp1Idx, 1.0 - fluidState.moleFraction(phase1Idx,comp0Idx));


        // set the additional quantities in the fluid state
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            paramCache.updateComposition(fluidState, phaseIdx);

            fluidState.setDensity(phaseIdx, FluidSystem::density(fluidState, paramCache, phaseIdx));
            fluidState.setMolarDensity(phaseIdx, FluidSystem::molarDensity(fluidState, paramCache, phaseIdx));
        }
    }
};

} // end namespace Dumux

#endif
