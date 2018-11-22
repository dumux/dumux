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
 * \brief Determines the pressures and saturations of all fluid phases
 *        given the total mass of all components.
 */
#ifndef DUMUX_COMPOSITIONAL_FLASH_HH
#define DUMUX_COMPOSITIONAL_FLASH_HH

#include <dune/common/fvector.hh>

#include <dumux/material/fluidstates/pseudo1p2c.hh>

namespace Dumux {

/*!
 * \ingroup ConstraintSolver
 * \brief Flash calculation routines for compositional sequential models
 *
 *        Routines for isothermal and isobaric 2p2c and 1p2c flash.
 */
template <class Scalar, class FluidSystem>
class CompositionalFlash
{
    using FluidState1p2c = PseudoOnePTwoCFluidState<Scalar, FluidSystem>;

    enum {  numPhases = FluidSystem::numPhases,
            numComponents = FluidSystem::numComponents
        };

    enum{
        phase0Idx = FluidSystem::phase0Idx,
        phase1Idx = FluidSystem::phase1Idx,
        comp0Idx = FluidSystem::comp0Idx,
        comp1Idx = FluidSystem::comp1Idx
    };

public:
    using ComponentVector = Dune::FieldVector<Scalar, numComponents>;
    using PhaseVector = Dune::FieldVector<Scalar, numPhases>;
/*!
 * \name Concentration flash for a given feed fraction
 */
//@{
    /*!
     * \brief 2p2c Flash for constant p & t if concentrations (feed mass fraction) is given.
     *
     * Routine goes as follows:
     * - determination of the equilibrium constants from the fluid system
     * - determination of maximum solubilities (mole fractions) according to phase pressures
     * - comparison with phase mass fraction Nu to determine phase presence => actual mole and mass fractions
     * - round off fluid properties
     * \param fluidState The sequential fluid State
     * \param Z0 Feed mass fraction: Mass of first component per total mass \f$\mathrm{[-]}\f$
     * \param phasePressure Vector holding the pressure \f$\mathrm{[Pa]}\f$
     * \param porosity Porosity \f$\mathrm{[-]}\f$
     * \param temperature Temperature \f$\mathrm{[K]}\f$
     */
    template<class FluidState>
    static void concentrationFlash2p2c(FluidState &fluidState,
                             const Scalar &Z0,
                             const PhaseVector &phasePressure,
                             const Scalar &porosity,
                             const Scalar &temperature)
    {
        // set the temperature, pressure
        fluidState.setTemperature(temperature);
        fluidState.setPressure(phase0Idx, phasePressure[phase0Idx]);
        fluidState.setPressure(phase1Idx, phasePressure[phase1Idx]);

        // mole equilibrium ratios k for in case first phase is reference phase
        Scalar k00 = FluidSystem::fugacityCoefficient(fluidState, phase0Idx, comp0Idx) * fluidState.pressure(phase0Idx)
                / (FluidSystem::fugacityCoefficient(fluidState, phase1Idx, comp0Idx) * fluidState.pressure(phase1Idx));
        Scalar k01 = FluidSystem::fugacityCoefficient(fluidState, phase0Idx, comp1Idx) * fluidState.pressure(phase0Idx)
                / (FluidSystem::fugacityCoefficient(fluidState, phase1Idx, comp1Idx) * fluidState.pressure(phase1Idx));

        // get mole fraction from equilibrium constants
        fluidState.setMoleFraction(phase0Idx,comp0Idx, ((1. - k01) / (k00 - k01)));
        fluidState.setMoleFraction(phase1Idx,comp0Idx, (fluidState.moleFraction(phase0Idx,comp0Idx) * k00));

        // transform mole to mass fractions
        fluidState.setMassFraction(phase0Idx, comp0Idx,
                (fluidState.moleFraction(phase0Idx, comp0Idx) * FluidSystem::molarMass(comp0Idx)
                / ( fluidState.moleFraction(phase0Idx, comp0Idx) * FluidSystem::molarMass(comp0Idx)
                    + (1.-fluidState.moleFraction(phase0Idx, comp0Idx)) * FluidSystem::molarMass(comp1Idx) )));
        fluidState.setMassFraction(phase1Idx, comp0Idx,
                (fluidState.moleFraction(phase1Idx, comp0Idx) * FluidSystem::molarMass(comp0Idx)
                / ( fluidState.moleFraction(phase1Idx, comp0Idx) * FluidSystem::molarMass(comp0Idx)
                    + (1.-fluidState.moleFraction(phase1Idx, comp0Idx)) * FluidSystem::molarMass(comp1Idx) )));

        // phase mass fraction Nu (ratio of phase mass to total phase mass) of first phase
        Scalar Nu0 = 1. + ((Z0 * (k00 - 1.)) + ((1. - Z0) * (k01 - 1.)))/((k01 - 1.) * (k00 - 1.));

        // check phase presence
        if (Nu0 > 0. && Nu0 < 1.) // two phases present
            fluidState.setNu(phase0Idx, Nu0);
        else if (Nu0 < 0.) // only second phase present
        {
            fluidState.setNu(phase0Idx, 0.);
            fluidState.setMassFraction(phase1Idx,comp0Idx, Z0); // assign complete mass dissolved into second phase

            // transform and store as mole fraction
            fluidState.setMoleFraction(phase1Idx,comp0Idx,
                                       (Z0 / FluidSystem::molarMass(comp0Idx)) / ((Z0 / FluidSystem::molarMass(comp0Idx)) + ((1- Z0) / FluidSystem::molarMass(comp1Idx))));
        }
        else // only first phase present
        {
            fluidState.setNu(phase0Idx, 1.); // no second phase
            fluidState.setMassFraction(phase0Idx, comp0Idx, Z0); // assign complete mass dissolved into first phase

            // transform and store as mole fraction
            fluidState.setMoleFraction(phase0Idx, comp0Idx,
                                       (Z0 / FluidSystem::molarMass(comp0Idx)) / ((Z0 / FluidSystem::molarMass(comp0Idx)) + ((1- Z0) / FluidSystem::molarMass(comp1Idx))));
        }

        // complete array of mass fractions
        fluidState.setMassFraction(phase0Idx, comp1Idx, 1. - fluidState.massFraction(phase0Idx,comp0Idx));
        fluidState.setMassFraction(phase1Idx, comp1Idx, 1. - fluidState.massFraction(phase1Idx,comp0Idx));
        // complete array of mole fractions
        fluidState.setMoleFraction(phase0Idx, comp1Idx, 1. - fluidState.moleFraction(phase0Idx,comp0Idx));
        fluidState.setMoleFraction(phase1Idx, comp1Idx, 1. - fluidState.moleFraction(phase1Idx,comp0Idx));

        // complete phase mass fractions
        fluidState.setNu(phase1Idx, 1. - fluidState.phaseMassFraction(phase0Idx));

        // get densities with correct composition
        fluidState.setDensity(phase0Idx, FluidSystem::density(fluidState, phase0Idx));
        fluidState.setDensity(phase1Idx, FluidSystem::density(fluidState, phase1Idx));
        fluidState.setMolarDensity(phase0Idx, FluidSystem::molarDensity(fluidState, phase0Idx));
        fluidState.setMolarDensity(phase1Idx, FluidSystem::molarDensity(fluidState, phase1Idx));

        Scalar sw = fluidState.phaseMassFraction(phase0Idx) / fluidState.density(phase0Idx);
        sw /= (fluidState.phaseMassFraction(phase0Idx)/fluidState.density(phase0Idx)
                    + fluidState.phaseMassFraction(phase1Idx)/fluidState.density(phase1Idx));
        fluidState.setSaturation(phase0Idx, sw);
    }

    /*!
     * \brief The simplest possible update routine for 1p2c "flash" calculations
     *
     * Routine goes as follows:
     * - Check if we are in single phase condition
     * - Assign total concentration to the present phase
     *
     * \param fluidState The sequential fluid state
     * \param Z0 Feed mass fraction: Mass of first component per total mass \f$\mathrm{[-]}\f$
     * \param phasePressure Vector holding the pressure \f$\mathrm{[Pa]}\f$
     * \param presentPhaseIdx Subdomain Index = Indication which phase is present
     * \param temperature Temperature \f$\mathrm{[K]}\f$
     */
    static void concentrationFlash1p2c(FluidState1p2c& fluidState, const Scalar& Z0,const Dune::FieldVector<Scalar,numPhases>
                                       phasePressure,const int presentPhaseIdx, const Scalar& temperature)
    {
        // set the temperature, pressure
        fluidState.setTemperature(temperature);
        fluidState.setPressure(phase0Idx, phasePressure[phase0Idx]);
        fluidState.setPressure(phase1Idx, phasePressure[phase1Idx]);

        fluidState.setPresentPhaseIdx(presentPhaseIdx);
        fluidState.setMassFraction(presentPhaseIdx,comp0Idx, Z0);

        // transform mass to mole fractions
        fluidState.setMoleFraction(presentPhaseIdx, comp0Idx, Z0 / FluidSystem::molarMass(comp0Idx)
                / (Z0 / FluidSystem::molarMass(comp0Idx) + (1. - Z0) / FluidSystem::molarMass(comp1Idx)));

        fluidState.setAverageMolarMass(presentPhaseIdx,
                fluidState.massFraction(presentPhaseIdx, comp0Idx) * FluidSystem::molarMass(comp0Idx)
                + fluidState.massFraction(presentPhaseIdx, comp1Idx) * FluidSystem::molarMass(comp1Idx));

        fluidState.setDensity(presentPhaseIdx, FluidSystem::density(fluidState, presentPhaseIdx));
        fluidState.setMolarDensity(presentPhaseIdx, FluidSystem::molarDensity(fluidState, presentPhaseIdx));
    }
//@}

/*!
 * \name Saturation flash for a given saturation (e.g. at boundary)
 */
//@{
    /*!
     * \brief A flash routine for 2p2c systems if the saturation instead of total concentration is known.
     *
     * Routine goes as follows:
     * - determination of the equilibrium constants from the fluid system
     * - determination of maximum solubilities (mole fractions) according to phase pressures
     * - round off fluid properties
     * \param fluidState The sequential fluid state
     * \param saturation Saturation of phase 1 \f$\mathrm{[-]}\f$
     * \param phasePressure Vector holding the pressure \f$\mathrm{[Pa]}\f$
     * \param porosity Porosity \f$\mathrm{[-]}\f$
     * \param temperature Temperature \f$\mathrm{[K]}\f$
     */
    template<class FluidState>
    static void saturationFlash2p2c(FluidState &fluidState,
            const Scalar &saturation,
            const PhaseVector &phasePressure,
            const Scalar &porosity,
            const Scalar &temperature)
    {
        // set the temperature, pressure
        fluidState.setTemperature(temperature);
        fluidState.setPressure(phase0Idx, phasePressure[phase0Idx]);
        fluidState.setPressure(phase1Idx, phasePressure[phase1Idx]);

        //mole equilibrium ratios K for in case wPhase is reference phase
        double k00 = FluidSystem::fugacityCoefficient(fluidState, phase0Idx, comp0Idx);
        double k01 = FluidSystem::fugacityCoefficient(fluidState, phase0Idx, comp1Idx);

        // get mole fraction from equilibrium constants
        fluidState.setMoleFraction(phase0Idx,comp0Idx, ((1. - k01) / (k00 - k01)));
        fluidState.setMoleFraction(phase1Idx,comp0Idx, (fluidState.moleFraction(phase0Idx,comp0Idx) * k00));

        // transform mole to mass fractions
        fluidState.setMassFraction(phase0Idx, comp0Idx,
                (fluidState.moleFraction(phase0Idx,comp0Idx) * FluidSystem::molarMass(comp0Idx)
                / ( fluidState.moleFraction(phase0Idx,comp0Idx) * FluidSystem::molarMass(comp0Idx)
                    + (1.-fluidState.moleFraction(phase0Idx,comp0Idx)) * FluidSystem::molarMass(comp1Idx) )));
        fluidState.setMassFraction(phase1Idx,comp0Idx,
                (fluidState.moleFraction(phase1Idx,comp0Idx) * FluidSystem::molarMass(comp0Idx)
                / ( fluidState.moleFraction(phase1Idx,comp0Idx) * FluidSystem::molarMass(comp0Idx)
                    + (1.-fluidState.moleFraction(phase1Idx,comp0Idx)) * FluidSystem::molarMass(comp1Idx) )));

        // complete array of mass fractions
        fluidState.setMassFraction(phase0Idx, comp1Idx, 1. - fluidState.massFraction(phase0Idx,comp0Idx));
        fluidState.setMassFraction(phase1Idx, comp1Idx, 1. - fluidState.massFraction(phase1Idx,comp0Idx));
        // complete array of mole fractions
        fluidState.setMoleFraction(phase0Idx, comp1Idx, 1. - fluidState.moleFraction(phase0Idx,comp0Idx));
        fluidState.setMoleFraction(phase1Idx, comp1Idx, 1. - fluidState.moleFraction(phase1Idx,comp0Idx));

        // get densities with correct composition
        fluidState.setDensity(phase0Idx, FluidSystem::density(fluidState, phase0Idx));
        fluidState.setDensity(phase1Idx, FluidSystem::density(fluidState, phase1Idx));
        fluidState.setMolarDensity(phase0Idx, FluidSystem::molarDensity(fluidState, phase0Idx));
        fluidState.setMolarDensity(phase1Idx, FluidSystem::molarDensity(fluidState, phase1Idx));

        // set saturation
        fluidState.setSaturation(phase0Idx, saturation);
    }
//@}
};

} // end namespace Dumux

#endif
