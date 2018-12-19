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
 * \ingroup TwoPOneCModel
 * \copydoc Dumux::TwoPOneCPrimaryVariableSwitch
 */

#ifndef DUMUX_2P1C_PRIMARY_VARIABLE_SWITCH_HH
#define DUMUX_2P1C_PRIMARY_VARIABLE_SWITCH_HH

#include <iostream>

#include <dumux/porousmediumflow/2p/formulation.hh>
#include <dumux/porousmediumflow/compositional/primaryvariableswitch.hh>

namespace Dumux {

/*!
 * \ingroup TwoPOneCModel
 * \brief The primary variable switch for the two-phase one-component model
 */
class TwoPOneCPrimaryVariableSwitch
: public PrimaryVariableSwitch<TwoPOneCPrimaryVariableSwitch>
{
    using ParentType = PrimaryVariableSwitch<TwoPOneCPrimaryVariableSwitch>;
    friend ParentType;

public:
    using ParentType::ParentType;

protected:

    /*!
     * \brief Performs variable switch at a degree of freedom location.
     *
     * \param priVars The primary variables at the given degree of freedom (dof) location.
     * \param volVars The volume variables.
     * \param dofIdxGlobal The respective dof index.
     * \param globalPos The global position of the dof.
     */
    template<class VolumeVariables, class GlobalPosition>
    bool update_(typename VolumeVariables::PrimaryVariables& priVars,
                 const VolumeVariables& volVars,
                 std::size_t dofIdxGlobal,
                 const GlobalPosition& globalPos)
    {
        using Scalar = typename VolumeVariables::PrimaryVariables::value_type;
        using FluidSystem = typename VolumeVariables::FluidSystem;
        using Indices = typename VolumeVariables::Indices;

        static constexpr auto formulation = VolumeVariables::priVarFormulation();
        static_assert( (formulation == TwoPFormulation::p0s1 || formulation == TwoPFormulation::p1s0),
                        "Chosen TwoPFormulation not supported!");

        // evaluate primary variable switch
        bool wouldSwitch = false;
        int phasePresence =  priVars.state();
        int newPhasePresence = phasePresence;

        // check if a primary var switch is necessary
        if (phasePresence == Indices::twoPhases)
        {
            Scalar Smin = 0.0;
            if (this->wasSwitched_[dofIdxGlobal])
                Smin = -0.01;

            if (volVars.saturation(FluidSystem::gasPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // gas phase disappears
                if (this->verbosity() > 1)
                    std::cout << "Gas phase (" << FluidSystem::phaseName(FluidSystem::gasPhaseIdx)
                              << ") disappears at dof " << dofIdxGlobal
                              << ", coordinates: " << globalPos
                              << ", S_" << FluidSystem::phaseName(FluidSystem::gasPhaseIdx) << ": "
                              << volVars.saturation(FluidSystem::gasPhaseIdx)
                              << std::endl;
                newPhasePresence = Indices::liquidPhaseOnly;

                priVars[Indices::switchIdx] = volVars.fluidState().temperature();
            }
            else if (volVars.saturation(FluidSystem::liquidPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // water phase disappears
                if (this->verbosity() > 1)
                    std::cout << "Liquid phase (" << FluidSystem::phaseName(FluidSystem::liquidPhaseIdx)
                              << ") disappears at dof " << dofIdxGlobal
                              << ", coordinates: " << globalPos
                              << ", S_" << FluidSystem::phaseName(FluidSystem::liquidPhaseIdx) << ": "
                              << volVars.saturation(FluidSystem::liquidPhaseIdx)
                              << std::endl;
                newPhasePresence = Indices::gasPhaseOnly;

                priVars[Indices::switchIdx] = volVars.fluidState().temperature();
            }

        }
        else if (phasePresence == Indices::liquidPhaseOnly)
        {
            const Scalar temp = volVars.fluidState().temperature();
            const Scalar tempVap = volVars.vaporTemperature();

            // if the the temperature would be larger than
            // the vapor temperature at the given pressure, gas phase appears
            if (temp >= tempVap)
            {
                wouldSwitch = true;
                // gas phase appears
                if (this->verbosity() > 1)
                    std::cout << "Gas phase (" << FluidSystem::phaseName(FluidSystem::gasPhaseIdx)
                              << ") appears at dof " << dofIdxGlobal
                              << ", coordinates: " << globalPos
                              << std::endl;
                newPhasePresence = Indices::twoPhases;
                if (formulation == TwoPFormulation::p1s0)
                    priVars[Indices::switchIdx] = 0.9999; // liquid phase saturation
                else
                    priVars[Indices::switchIdx] = 0.0001;
            }
        }
        else if (phasePresence == Indices::gasPhaseOnly)
        {
            const Scalar temp = volVars.fluidState().temperature();
            const Scalar tempVap = volVars.vaporTemperature();

            if (temp < tempVap)
            {
                wouldSwitch = true;
                // liquid phase appears
                if (this->verbosity() > 1)
                    std::cout << "Liquid phase (" << FluidSystem::phaseName(FluidSystem::liquidPhaseIdx) << ") appears at dof " << dofIdxGlobal
                              << ", coordinates: " << globalPos  << std::endl;

               newPhasePresence = Indices::twoPhases;
               if (formulation == TwoPFormulation::p1s0)
                   priVars[Indices::switchIdx] = 0.0001;
               else
                   priVars[Indices::switchIdx] = 0.9999;
            }
    }
        priVars.setState(newPhasePresence);
        this->wasSwitched_[dofIdxGlobal] = wouldSwitch;
        return phasePresence != newPhasePresence;
    }
};

} // end namespace Dumux

#endif
