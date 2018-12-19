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
 * \ingroup RichardsModel
 * \brief The primary variable switch for the extended Richards model.
 */

#ifndef DUMUX_RICHARDS_PRIMARY_VARIABLE_SWITCH_HH
#define DUMUX_RICHARDS_PRIMARY_VARIABLE_SWITCH_HH

#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>
#include <dumux/material/constants.hh>
#include <dumux/porousmediumflow/compositional/primaryvariableswitch.hh>

namespace Dumux {

/*!
 * \ingroup RichardsModel
 * \brief The primary variable switch controlling the phase presence state variable.
 */
class ExtendedRichardsPrimaryVariableSwitch
: public PrimaryVariableSwitch<ExtendedRichardsPrimaryVariableSwitch>
{
    using ParentType = PrimaryVariableSwitch<ExtendedRichardsPrimaryVariableSwitch>;
    friend ParentType;

public:
    using ParentType::ParentType;

protected:

    // perform variable switch at a degree of freedom location
    template<class VolumeVariables, class GlobalPosition>
    bool update_(typename VolumeVariables::PrimaryVariables& priVars,
                 const VolumeVariables& volVars,
                 std::size_t dofIdxGlobal,
                 const GlobalPosition& globalPos)
    {
        using Scalar = typename VolumeVariables::PrimaryVariables::value_type;
        using Indices = typename VolumeVariables::Indices;
        using FluidSystem = typename VolumeVariables::FluidSystem;

        static const bool usePriVarSwitch = getParam<bool>("Problem.UsePrimaryVariableSwitch");
        if (!usePriVarSwitch)
            return false;

        if (!VolumeVariables::enableWaterDiffusionInAir())
            DUNE_THROW(Dune::InvalidStateException, "The Richards primary variable switch only works with water diffusion in air enabled!");

        static constexpr int liquidCompIdx = FluidSystem::liquidPhaseIdx;

        // evaluate primary variable switch
        bool wouldSwitch = false;
        int phasePresence = priVars.state();
        int newPhasePresence = phasePresence;

        // check if a primary var switch is necessary
        if (phasePresence == Indices::gasPhaseOnly)
        {
            // if the mole fraction of water is larger than the one
            // predicted by a liquid-vapor equilibrium
            Scalar xnw = volVars.moleFraction(FluidSystem::gasPhaseIdx, liquidCompIdx);
            Scalar xnwPredicted = FluidSystem::H2O::vaporPressure(volVars.temperature())
                                  / volVars.pressure(FluidSystem::gasPhaseIdx);

            Scalar xwMax = 1.0;
            if (xnw / xnwPredicted > xwMax)
                wouldSwitch = true;
            if (this->wasSwitched_[dofIdxGlobal])
                xwMax *= 1.01;

            // if the ratio of predicted mole fraction to current mole fraction is larger than
            // 100%, wetting phase appears
            if (xnw / xnwPredicted > xwMax)
            {
                // wetting phase appears
                if (this->verbosity() > 1)
                    std::cout << "Liquid phase appears at dof " << dofIdxGlobal
                              << ", coordinates: " << globalPos << ", xnw / xnwPredicted * 100: "
                              << xnw / xnwPredicted * 100 << "%"
                              << ", at x_n^w: " << priVars[Indices::switchIdx] << std::endl;
                newPhasePresence = Indices::bothPhases;
                priVars[Indices::switchIdx] = 0.0;
            }
        }
        else if (phasePresence == Indices::bothPhases)
        {
            Scalar Smin = 0.0;
            if (this->wasSwitched_[dofIdxGlobal])
                Smin = -0.01;

            if (volVars.saturation(FluidSystem::liquidPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // wetting phase disappears
                newPhasePresence = Indices::gasPhaseOnly;
                priVars[Indices::switchIdx] = volVars.moleFraction(FluidSystem::gasPhaseIdx, liquidCompIdx);

                if (this->verbosity() > 1)
                    std::cout << "Liquid phase disappears at dof " << dofIdxGlobal
                              << ", coordinates: " << globalPos << ", sw: "
                              << volVars.saturation(FluidSystem::liquidPhaseIdx)
                              << ", x_n^w: " << priVars[Indices::switchIdx] << std::endl;
            }
        }
        else if (phasePresence == Indices::liquidPhaseOnly)
        {
            DUNE_THROW(Dune::NotImplemented, "Water phase only phase presence!");
        }

        priVars.setState(newPhasePresence);
        this->wasSwitched_[dofIdxGlobal] = wouldSwitch;
        return phasePresence != newPhasePresence;
    }
};

} // end namespace Dumux

#endif
