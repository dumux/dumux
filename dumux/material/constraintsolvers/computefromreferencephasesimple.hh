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
 * \brief Computes the equilibrium composition of a generic fluid state if a
 *        reference phase has been specified.
 */
#ifndef DUMUX_COMPUTE_FROM_REFERENCE_PHASE_SIMPLE_HH
#define DUMUX_COMPUTE_FROM_REFERENCE_PHASE_SIMPLE_HH

#include <dumux/common/exceptions.hh>
#include <dumux/common/valgrind.hh>

namespace Dumux {

/*!
 * \ingroup ConstraintSolvers
 * \brief Computes the equilibrium composition of a generic fluid state if a
 *        reference phase has been specified.
 *
 * This makes it possible to specify just one phase and let the
 * composition of the remaining ones be calculated by the constraint solver.
 * This explicit compositional flash can only handle ideal mixtures
 * and only works for two-phase, two-component and three-phase,
 * three-component systems.
 * It assumes the following quantities to be set:
 *
 * - composition (mole/mass fractions) of the *reference* phase
 *   \f$x^\kappa_\beta\f$, \f$X^\kappa_\beta\f$
 * - temperature of the *reference* phase \f$T_\beta\f$
 * - saturations of *all* phases (i.e., all are present)
 * - pressures of *all* phases \f$p_\alpha\f$, \f$p_\beta\f$
 *
 * After calling the solve() method the following quantities are
 * calculated:
 *
 * - fugacity coefficients of *all* components in *all* phases
 *   \f$\Phi^\kappa_\alpha\f$, \f$\Phi^\kappa_\beta\f$
 * - composition in mole fractions of *all* components in *all* phases \f$x^\kappa_\alpha\f$
 * - density and molar density of *all* phases
 *   \f$\rho_\alpha\f$, \f$\rho_\beta\f$, \f$\rho_{mol, \alpha}\f$, \f$\rho_{mol, \beta}\f$
 *
 * Setting those quantities also guarantees by design of the fluid state
 * that the following quantities are set:
 *
 * - temperature of *all* phases \f$T_\alpha\f$, \f$T_\beta\f$
 * - molar volume of *all* phases
 *   \f$V_{mol, \alpha}\f$, \f$V_{mol, \beta}\f$
 * - mass fractions and molarities of *all* components in *all* phases
 *   \f$X^\kappa_\alpha\f$, \f$X^\kappa_\beta\f$,
 *   \f$c^\kappa_\alpha\f$, \f$c^\kappa_\beta\f$
 * - mean molar masses of *all* phases \f$M_\alpha\f$, \f$M_\beta\f$
 */
template <class Scalar, class FluidSystem>
class ComputeFromReferencePhaseSimple
{
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };

    // component and phase indices
    enum
    {
        comp0Idx = FluidSystem::comp0Idx,
        comp1Idx = FluidSystem::comp1Idx,
        phase0Idx = FluidSystem::phase0Idx,
        phase1Idx = FluidSystem::phase1Idx
    };

public:
    /*!
     * \brief Computes all quantities of a generic fluid state if a
     *        reference phase has been specified.
     *
     * \param fluidState Thermodynamic state of the fluids
     * \param paramCache  Container for cache parameters
     * \param refPhaseIdx The phase index of the reference phase
     */
    template <class FluidState, class ParameterCache>
    static void solve(FluidState &fluidState,
                      ParameterCache &paramCache,
                      int refPhaseIdx)
    {
        // restrictions of the simple compositional flash method
#ifndef NDEBUG
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            assert(FluidSystem::isIdealMixture(phaseIdx));
        }
#endif
        static_assert((numPhases == 2 && numComponents == 2) || (numPhases == 3 && numComponents == 3),
                      "The simple compositional flash can only handle two-phase, two-component or three-phase, three-component system");

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

            const Scalar xref0 = fluidState.moleFraction(refPhaseIdx, comp0Idx);
            const Scalar xref1 = fluidState.moleFraction(refPhaseIdx, comp1Idx);

            // get phase pressures
            const Scalar pl = fluidState.pressure(liquidPhaseIdx);
            const Scalar pg = fluidState.pressure(gasPhaseIdx);

            // note that the non-reference phase is not actually existing!
            // Thus, this is used as phase switch criterion.
            if(refPhaseIdx == liquidPhaseIdx)
            {
                // Dalton's law is used to calculate equilibrium mole fractions in the gas phase:
                // \f$ x_g^{\kappa} = p_g^{\kappa} / p_g \f$.
                // The partial pressure of the liquid component in the gas phase is calculated using Henry's law:
                // \f$ x_l^l = p_g^l / H^l \f$.
                // The partial pressure of the gas component in the gas phase is calculated using Raoult's law:
                // \f$ p^g_g = p_{vap}^g x_l^g \f$.

                const Scalar xg0 = xref0*( FluidSystem::fugacityCoefficient(fluidState, liquidPhaseIdx, comp0Idx)*pl )/pg;
                const Scalar xg1 = xref1*( FluidSystem::fugacityCoefficient(fluidState, liquidPhaseIdx, comp1Idx)*pl )/pg;

                fluidState.setMoleFraction(gasPhaseIdx, comp0Idx, xg0);
                fluidState.setMoleFraction(gasPhaseIdx, comp1Idx, xg1);
            }
            else if(refPhaseIdx == phase1Idx)
            {
                // Henry's law is used to calculate equilibrium mole fraction of the liquid component in the liquid phase.
                // Raoult's law is used to calculate equilibrium mole fraction of the gas component in the liquid phase.
                // Dalton's law is used to calculate the partial pressures in the gas phase.

                const Scalar xl0 = xref0*pg/( FluidSystem::fugacityCoefficient(fluidState, liquidPhaseIdx, comp0Idx)*pl );
                const Scalar xl1 = xref1*pg/( FluidSystem::fugacityCoefficient(fluidState, liquidPhaseIdx, comp1Idx)*pl );

                fluidState.setMoleFraction(liquidPhaseIdx, comp0Idx, xl0);
                fluidState.setMoleFraction(liquidPhaseIdx, comp1Idx, xl1);
            }
        }
        else if(numPhases == 3)
        {
            DUNE_THROW(Dune::NotImplemented, "3p3c not yet implemented");
        }

        // set additional secondary variables in the fluid state for all phases
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
