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
 * \ingroup CO2Model
 * \brief The primary variable switch for the 2p2c-CO2 model
 */

#ifndef DUMUX_2P2C_CO2_PRIMARY_VARIABLE_SWITCH_HH
#define DUMUX_2P2C_CO2_PRIMARY_VARIABLE_SWITCH_HH

#include <iostream>

#include <dumux/porousmediumflow/compositional/primaryvariableswitch.hh>
#include <dumux/porousmediumflow/2p/formulation.hh>

namespace Dumux
{
/*!
 * \ingroup CO2Model
 * \brief The primary variable switch for the 2p2c-CO2 model controlling the phase presence state variable.
 *
 * The phase switch occurs when the equilibrium concentration
 * of a component in a phase is exceeded, instead of the sum of the components in the virtual phase
 * (the phase which is not present) being greater that unity as done in the 2p2c model.
 */
class TwoPTwoCCO2PrimaryVariableSwitch
: public PrimaryVariableSwitch< TwoPTwoCCO2PrimaryVariableSwitch >
{
    using ParentType = PrimaryVariableSwitch< TwoPTwoCCO2PrimaryVariableSwitch >;
    friend ParentType;

public:
    using ParentType::ParentType;

protected:
    // perform variable switch at a degree of freedom location
    template<class VolumeVariables, class IndexType, class GlobalPosition>
    bool update_(typename VolumeVariables::PrimaryVariables& priVars,
                 const VolumeVariables& volVars,
                 IndexType dofIdxGlobal,
                 const GlobalPosition& globalPos)
    {
        using Scalar = typename VolumeVariables::PrimaryVariables::value_type;

        using FluidSystem = typename VolumeVariables::FluidSystem;
        static constexpr int phase0Idx = FluidSystem::phase0Idx;
        static constexpr int phase1Idx = FluidSystem::phase1Idx;
        static constexpr int comp0Idx = FluidSystem::comp0Idx;
        static constexpr int comp1Idx = FluidSystem::comp1Idx;

        static constexpr bool useMoles = VolumeVariables::useMoles();
        static constexpr auto formulation = VolumeVariables::priVarFormulation();
        static_assert( (formulation == TwoPFormulation::p0s1 || formulation == TwoPFormulation::p1s0),
                        "Chosen TwoPFormulation not supported!");

        using Indices = typename VolumeVariables::Indices;
        static constexpr int switchIdx = Indices::switchIdx;

        // evaluate primary variable switch
        bool wouldSwitch = false;
        int phasePresence = priVars.state();
        int newPhasePresence = phasePresence;

        // the param cache to evaluate the equilibrium mole fraction
        typename FluidSystem::ParameterCache paramCache;

        // check if a primary var switch is necessary
        if (phasePresence == Indices::secondPhaseOnly)
        {
            // calculate wetting component mole fraction in the second phase
            Scalar xnw = volVars.moleFraction(phase1Idx, comp0Idx);
            Scalar xnwMax = FluidSystem::equilibriumMoleFraction(volVars.fluidState(), paramCache, phase1Idx);

            // if it is larger than the equilibirum mole fraction switch
            if(xnw > xnwMax)
                wouldSwitch = true;

            if (this->wasSwitched_[dofIdxGlobal])
                xnwMax *= 1.02;

            // if it is larger than the equilibirum mole fraction switch: first phase appears
            if (xnw > xnwMax)
            {
                // wetting phase appears
                if (this->verbosity() > 1)
                    std::cout << "First phase (" << FluidSystem::phaseName(phase0Idx) << ") appears at dof " << dofIdxGlobal
                              << ", coordinates: " << globalPos
                              << ", x^" << FluidSystem::componentName(comp0Idx) << "_" << FluidSystem::phaseName(phase1Idx) << " > x_equilibrium: "
                              << xnw << " > " << xnwMax << std::endl;
                newPhasePresence = Indices::bothPhases;
                if (formulation == TwoPFormulation::p1s0)
                    priVars[switchIdx] = 0.0;
                else
                    priVars[switchIdx] = 1.0;
            }
        }
        else if (phasePresence == Indices::firstPhaseOnly)
        {
            // calculate second component mole fraction in the wetting phase
            Scalar xwn = volVars.moleFraction(phase0Idx, comp1Idx);
            Scalar xwnMax = FluidSystem::equilibriumMoleFraction(volVars.fluidState(), paramCache, phase0Idx);

            // if it is larger than the equilibirum mole fraction switch
            if(xwn > xwnMax)
                wouldSwitch = true;

            if (this->wasSwitched_[dofIdxGlobal])
                xwnMax *= 1.02;

            // if it is larger than the equilibirum mole fraction switch second phase appears
            if(xwn > xwnMax)
            {
                // Second phase appears
                if (this->verbosity() > 1)
                    std::cout << "Second phase (" << FluidSystem::phaseName(phase1Idx) << ") appears at dof " << dofIdxGlobal
                              << ", coordinates: " << globalPos
                              << ", x^" << FluidSystem::componentName(comp1Idx) << "_" << FluidSystem::phaseName(phase0Idx) << " > x_equilibrium: "
                              << xwn << " > " << xwnMax << std::endl;
                newPhasePresence = Indices::bothPhases;
                if (formulation == TwoPFormulation::p1s0)
                    priVars[switchIdx] = 0.999;
                else
                    priVars[switchIdx] = 0.001;
            }
        }
        // TODO: this is the same as for the 2p2c model maybe factor out
        else if (phasePresence == Indices::bothPhases)
        {
            Scalar Smin = 0.0;
            if (this->wasSwitched_[dofIdxGlobal])
                Smin = -0.01;

            if (volVars.saturation(phase1Idx) <= Smin)
            {
                wouldSwitch = true;
                // nonwetting phase disappears
                if (this->verbosity() > 1)
                    std::cout << "Second phase (" << FluidSystem::phaseName(phase1Idx) << ") disappears at dof " << dofIdxGlobal
                              << ", coordinates: " << globalPos
                              << ", S_" << FluidSystem::phaseName(phase1Idx) << ": " << volVars.saturation(phase1Idx)
                              << std::endl;
                newPhasePresence = Indices::firstPhaseOnly;

                if(useMoles) // mole-fraction formulation
                    priVars[switchIdx] = volVars.moleFraction(phase0Idx, comp1Idx);
                else // mass-fraction formulation
                    priVars[switchIdx] = volVars.massFraction(phase0Idx, comp1Idx);
            }
            else if (volVars.saturation(phase0Idx) <= Smin)
            {
                wouldSwitch = true;
                // wetting phase disappears
                if (this->verbosity() > 1)
                    std::cout << "First phase (" << FluidSystem::phaseName(phase0Idx) << ") disappears at dof " << dofIdxGlobal
                              << ", coordinates: " << globalPos
                              << ", S_" << FluidSystem::phaseName(phase0Idx) << ": " << volVars.saturation(phase0Idx)
                              << std::endl;
                newPhasePresence = Indices::secondPhaseOnly;

                if(useMoles) // mole-fraction formulation
                    priVars[switchIdx] = volVars.moleFraction(phase1Idx, comp0Idx);
                else // mass-fraction formulation
                    priVars[switchIdx] = volVars.massFraction(phase1Idx, comp0Idx);
            }
        }

        priVars.setState(newPhasePresence);
        this->wasSwitched_[dofIdxGlobal] = wouldSwitch;
        return phasePresence != newPhasePresence;
    }
};

} // end namespace Dumux

#endif
