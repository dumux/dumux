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
 * \brief Computes the equilibrium composition of all phases
 */
#ifndef DUMUX_MISCIBLE_MULTIPHASE_COMPOSITION_SIMPLE_HH
#define DUMUX_MISCIBLE_MULTIPHASE_COMPOSITION_SIMPLE_HH

#include <dumux/common/exceptions.hh>
#include <dumux/common/valgrind.hh>

namespace Dumux {
/*!
 * \ingroup ConstraintSolvers
 * \brief Computes the equilibrium composition of all phases
 *
 * This 'simple', explicit compositional flash can only handle ideal mixtures
 * and only works for two-phase, two-component and three-phase,
 * three-component systems.
 * It assumes the following quantities to be set:
 *
 * - temperatures of *all* phases \f$T_\alpha\f$
 * - pressures of *all* phases \f$p_\alpha\f$
 * - saturations of *all* phases (i.e., all are present)
 *
 * It also assumes that the mole/mass fractions of all phases sum up
 * to 1. Henry's law is used to calculate the composition of the liquid phase:
 * \f$ x_l^{\kappa} = p_g^{\kappa} / H^g \f$.
 * Dalton's law is used to calculate the composition of the gas phase:
 * \f$ x_g^{\kappa} = p_g^{\kappa} / p_g \f$.
 * Here, the unknown partial pressures
 * are approximated using the vapor pressure of the liquid component
 * \f$ p_g^l \approx p_{vap}^l\f$ (Raoult's law) and substituting
 * \f$ p_g^g = p_g - p_g^l\f$.
 *
 * After calling the solve() method, the following quantities
 * are calculated:
 *
 * - 'fugacity coefficients' of *all* components in *all* phases
 *   (Note that this is done for the sake of harmonisation
 *   with the other flashs. This flash does not use fugacity coefficients.)
 *   \f$\Phi^\kappa_\alpha\f$
 * - composition in mole fractions of *all* components in *all* phases \f$x^\kappa_\alpha\f$
 * - density and molar density of *all* phases
 *   \f$\rho_\alpha\f$, \f$\rho_{mol, \alpha}\f$
 *
 * Setting those quantities also guarantees by design of the fluid state
 * that the following quantities are set:
 *
 * - molar volume of *all* phases \f$V_{mol, \alpha}\f$
 * - mean molar masses of *all* phases \f$M_\alpha\f$
 * - mass fractions and molarities of *all* components in *all* phases
 *   \f$X^\kappa_\alpha\f$, \f$c^\kappa_\alpha\f$
 */
template <class Scalar, class FluidSystem>
class MiscibleMultiPhaseCompositionSimple
{
    static constexpr int numPhases = FluidSystem::numPhases;
    static constexpr int numComponents = FluidSystem::numComponents;

    // phase indices
    enum
    {
        phase0Idx = FluidSystem::phase0Idx,
        phase1Idx = FluidSystem::phase1Idx
    };

public:
    /*!
     * \brief @copybrief Dumux::MiscibleMultiPhaseCompositionSimple
     *
     * \param fluidState A container with the current (physical) state of the fluid
     * \param paramCache A container for iterative calculation of fluid composition
     * \param knownPhaseIdx The index of the phase with known properties
     */
    template <class FluidState, class ParameterCache>
    static void solve(FluidState &fluidState,
                      ParameterCache &paramCache,
                      int knownPhaseIdx = 0)
    {
        // restrictions of the simple compositional flash method
#ifndef NDEBUG
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            assert(FluidSystem::isIdealMixture(phaseIdx));
        }
#endif
        static_assert((numPhases == 2 && numComponents == 2) || (numPhases == 3 && numComponents == 3),
                      "The simple compositional flash can only handle two-phase, two-component or three-phase, three-component systems");

        // set all fugacity coefficients
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx)
        {
            for (int compIdx = 0; compIdx < numComponents; ++ compIdx)
            {
                Scalar phi = FluidSystem::fugacityCoefficient(fluidState, paramCache, phaseIdx, compIdx);
                fluidState.setFugacityCoefficient(phaseIdx, compIdx, phi);
            }
        }

        if(numPhases == 2)
        {
            // determine indices
            int liquidPhaseIdx;
            int gasPhaseIdx;
            if (FluidSystem::isGas(phase1Idx))
            {
                liquidPhaseIdx = phase0Idx;
                gasPhaseIdx = phase1Idx;
            }
            else if (FluidSystem::isGas(phase0Idx))
            {
                liquidPhaseIdx = phase1Idx;
                gasPhaseIdx = phase0Idx;
            }
            else
            {
                DUNE_THROW(Dune::InvalidStateException, "The simple compositional flash assumes one gas and one liquid phase");
            }
            static const int liquidCompIdx = FluidSystem::getMainComponent(liquidPhaseIdx);
            static const int gasCompIdx = FluidSystem::getMainComponent(gasPhaseIdx);

            // get phase pressures
            const Scalar pl = fluidState.pressure(liquidPhaseIdx);
            const Scalar pg = fluidState.pressure(gasPhaseIdx);

            // use Daltons's law to determine mole fraction of the gas component in the gas phase:
            // \f$ x_g^g = p_g^g / p_g  = (p_g - p_g^l) / p_g \f$.
            // We assume that the partial pressure of the liquid component in the gas phase can be
            // approximated with the vapor pressure of the liquid component.
            // This neglects the influence of other components in the gas phase.
            const Scalar partPressLiquidCompInGas = FluidSystem::vaporPressure(fluidState, liquidCompIdx);

            // get the partial pressure of the gas component in the gas phase
            const Scalar partPressGasCompInGas = pg - partPressLiquidCompInGas;

            // calculate the mole fraction of the gas component in the gas phase
            Scalar xgg = partPressGasCompInGas / pg;
            xgg = std::min(1.0, std::max(0.0, xgg)); // regularization

            // calculate the mole fraction of the liquid component in the gas phase
            const Scalar xgl = 1.0 - xgg;

            // use Henry's law to determine mole fraction of the gas component in the liquid phase
            // \f$ x_l^g = p_g^g / H^g \f$.
            Scalar xlg = partPressGasCompInGas / FluidSystem::henry(fluidState, liquidPhaseIdx, gasCompIdx);
            xlg = std::min(1.0, std::max(0.0, xlg)); // regularization

            // calculate the mole fraction of the liquid component in the liquid phase:
            const Scalar xll = 1.0 - xlg;

            // set all mole fractions
            fluidState.setMoleFraction(liquidPhaseIdx, liquidCompIdx, xll);
            fluidState.setMoleFraction(liquidPhaseIdx, gasCompIdx, xlg);
            fluidState.setMoleFraction(gasPhaseIdx, liquidCompIdx, xgl);
            fluidState.setMoleFraction(gasPhaseIdx, gasCompIdx, xgg);
        }
        else if(numPhases == 3)
        {
            DUNE_THROW(Dune::NotImplemented, "3p3c not yet implemented");
        }

        // set additional secondary variables in the fluid state
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            paramCache.updateComposition(fluidState, phaseIdx);

            Scalar value = FluidSystem::density(fluidState, paramCache, phaseIdx);
            fluidState.setDensity(phaseIdx, value);

            value = FluidSystem::molarDensity(fluidState, paramCache, phaseIdx);
            fluidState.setMolarDensity(phaseIdx, value);
        }
    }
};

} // end namespace Dumux

#endif
