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
 * \ingroup TwoPNCModel
 * \brief The primary variable switch for the 2pnc model.
 */

#ifndef DUMUX_2PNC_PRIMARY_VARIABLE_SWITCH_HH
#define DUMUX_2PNC_PRIMARY_VARIABLE_SWITCH_HH

#include <iostream>

#include <dumux/porousmediumflow/compositional/primaryvariableswitch.hh>
#include <dumux/porousmediumflow/2p/formulation.hh>

namespace Dumux {

/*!
 * \ingroup TwoPNCModel
 * \brief The primary variable switch controlling the phase presence state variable.
 */
class TwoPNCPrimaryVariableSwitch
: public PrimaryVariableSwitch<TwoPNCPrimaryVariableSwitch>
{
    using ParentType = PrimaryVariableSwitch<TwoPNCPrimaryVariableSwitch>;
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

        static constexpr auto numComponents = VolumeVariables::numFluidComponents();
        static constexpr bool useMoles = VolumeVariables::useMoles();
        static_assert(useMoles || numComponents < 3, "!useMoles is only implemented for numComponents < 3.");
        static constexpr auto numMajorComponents = VolumeVariables::numFluidPhases();
        static constexpr auto formulation = VolumeVariables::priVarFormulation();
        static_assert( (formulation == TwoPFormulation::p0s1 || formulation == TwoPFormulation::p1s0),
                        "Chosen TwoPFormulation not supported!");

        using Indices = typename VolumeVariables::Indices;
        static constexpr int switchIdx = Indices::switchIdx;

        // evaluate primary variable switch
        bool wouldSwitch = false;
        int phasePresence = priVars.state();
        int newPhasePresence = phasePresence;

        //check if a primary variable switch is necessary
        if (phasePresence == Indices::bothPhases)
        {
            Scalar Smin = 0; //saturation threshold
            if (this->wasSwitched_[dofIdxGlobal])
                Smin = -0.01;

            // if saturation of first phase is smaller 0: switch
            if (volVars.saturation(phase0Idx) <= Smin)
            {
                wouldSwitch = true;
                // first phase has to disappear
                if (this->verbosity() > 1)
                    std::cout << "First phase (" << FluidSystem::phaseName(phase0Idx) << ")"
                              << " disappears at dof " << dofIdxGlobal
                              << ", coordinates: " << globalPos
                              << ", S_" << FluidSystem::phaseName(phase0Idx) << ": " << volVars.saturation(phase0Idx)
                              << std::endl;
                newPhasePresence = Indices::secondPhaseOnly;

                // switch not depending on formulation, switch "S0" to "x10"
                if(useMoles) // mole-fraction formulation
                    priVars[switchIdx] = volVars.moleFraction(phase1Idx, comp0Idx);
                else // mass-fraction formulation
                    priVars[switchIdx] = volVars.massFraction(phase1Idx, comp0Idx);

                // switch all secondary components to mole fraction in nonwetting phase
                if(useMoles) // mole-fraction formulation
                    for (int compIdx = numMajorComponents; compIdx < numComponents; ++compIdx)
                        priVars[compIdx] = volVars.moleFraction(phase1Idx, compIdx);
                else // mass-fraction formulation
                    for (int compIdx = numMajorComponents; compIdx < numComponents; ++compIdx)
                        priVars[compIdx] = volVars.massFraction(phase1Idx, compIdx);
            }

            // if saturation of second phase is smaller than 0: switch
            else if (volVars.saturation(phase1Idx) <= Smin)
            {
                wouldSwitch = true;
                // second phase has to disappear
                if (this->verbosity() > 1)
                    std::cout << "Second phase (" << FluidSystem::phaseName(phase1Idx) << ")"
                              << " disappears at dof " << dofIdxGlobal
                              << ", coordinates: " << globalPos
                              << ", S_" << FluidSystem::phaseName(phase1Idx) << ": " << volVars.saturation(phase1Idx)
                              << std::endl;
                newPhasePresence = Indices::firstPhaseOnly;

                // switch "S1" to "x01"
                if(useMoles) // mole-fraction formulation
                    priVars[switchIdx] = volVars.moleFraction(phase0Idx, comp1Idx);
                else // mass-fraction formulation
                    priVars[switchIdx] = volVars.massFraction(phase0Idx, comp1Idx);
            }
        }
        else if (phasePresence == Indices::secondPhaseOnly)
        {
            Scalar x0Max = 1;
            Scalar x0Sum = 0;
            // Calculate sum of mole fractions in the hypothetical first phase
            for (int compIdx = 0; compIdx < numComponents; compIdx++)
                x0Sum += volVars.moleFraction(phase0Idx, compIdx);

            if (x0Sum > x0Max)
                wouldSwitch = true;
            if (this->wasSwitched_[dofIdxGlobal])
                x0Max *= 1.02;

            // first phase appears if sum is larger than one
            if (x0Sum/*sum of mole fractions*/ > x0Max/*1*/)
            {
                if (this->verbosity() > 1)
                    std::cout << "Second phase (" << FluidSystem::phaseName(phase0Idx) << ")"
                              << " appears at dof " << dofIdxGlobal
                              << ", coordinates: " << globalPos
                              << ", sum x^i_" << FluidSystem::phaseName(phase0Idx) << ": " << x0Sum
                              << std::endl;
                newPhasePresence = Indices::bothPhases;

                // saturation of the first phase set to 0.0001 (if formulation TwoPFormulation::p1s0 and vice versa)
                if (formulation == TwoPFormulation::p1s0)
                    priVars[switchIdx] = 0.0001;
                else
                    priVars[switchIdx] = 0.9999;

                // switch all secondary components back to first component mole fraction
                for (int compIdx = numMajorComponents; compIdx < numComponents; ++compIdx)
                    priVars[compIdx] = volVars.moleFraction(phase0Idx,compIdx);
            }
        }
        else if (phasePresence == Indices::firstPhaseOnly)
        {
            Scalar x1Max = 1;
            Scalar x1Sum = 0;

            // Calculate sum of mole fractions in the hypothetical wetting phase
            for (int compIdx = 0; compIdx < numComponents; compIdx++)
                x1Sum += volVars.moleFraction(phase1Idx, compIdx);

            if (x1Sum > x1Max)
                wouldSwitch = true;
            if (this->wasSwitched_[dofIdxGlobal])
                x1Max *= 1.02;

            // wetting phase appears if sum is larger than one
            if (x1Sum > x1Max)
            {
                if (this->verbosity() > 1)
                    std::cout << "Second phase (" << FluidSystem::phaseName(phase1Idx) << ")"
                              << " appears at dof " << dofIdxGlobal
                              << ", coordinates: " << globalPos
                              << ", sum x^i_" << FluidSystem::phaseName(phase1Idx) << ": " << x1Sum
                              << std::endl;
                newPhasePresence = Indices::bothPhases;
                //saturation of the wetting phase set to 0.9999 (if formulation TwoPFormulation::pnsw and vice versa)
                if (formulation == TwoPFormulation::p1s0)
                    priVars[switchIdx] = 0.9999;
                else
                    priVars[switchIdx] = 0.0001;
            }
        }

        priVars.setState(newPhasePresence);
        this->wasSwitched_[dofIdxGlobal] = wouldSwitch;
        return phasePresence != newPhasePresence;
    }
};

} // end namespace Dumux

#endif
