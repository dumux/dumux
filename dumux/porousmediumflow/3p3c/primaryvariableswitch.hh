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
 * \ingroup ThreePThreeCModel
 * \brief The primary variable switch for the 3p3c model.
 */

#ifndef DUMUX_3P3C_PRIMARY_VARIABLE_SWITCH_HH
#define DUMUX_3P3C_PRIMARY_VARIABLE_SWITCH_HH

#include <iostream>

#include <dumux/porousmediumflow/compositional/primaryvariableswitch.hh>

namespace Dumux {

/*!
 * \ingroup ThreePThreeCModel
 * \brief The primary variable switch controlling the phase presence state variable.
 */
class ThreePThreeCPrimaryVariableSwitch
: public PrimaryVariableSwitch<ThreePThreeCPrimaryVariableSwitch>
{
    using ParentType = PrimaryVariableSwitch<ThreePThreeCPrimaryVariableSwitch>;
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
        using PrimaryVariables = typename VolumeVariables::PrimaryVariables;
        using Scalar = typename PrimaryVariables::value_type;
        using Indices = typename VolumeVariables::Indices;
        using FluidSystem = typename VolumeVariables::FluidSystem;

        // evaluate if the primary variable switch would switch
        bool wouldSwitch = false;
        auto phasePresence = priVars.state();
        auto newPhasePresence = phasePresence;

        // check if a primary variable switch is necessary
        if (phasePresence == Indices::threePhases)
        {
            Scalar Smin = 0;
            if (this->wasSwitched_[dofIdxGlobal])
                Smin = -0.01;

            if (volVars.saturation(FluidSystem::gPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // gas phase disappears
                if (this->verbosity() > 1)
                    std::cout << "Gas phase disappears at dof " << dofIdxGlobal
                              << ", coordinates: " << globalPos << ", sg: "
                              << volVars.saturation(FluidSystem::gPhaseIdx) << std::endl;
                newPhasePresence = Indices::wnPhaseOnly;

                priVars[Indices::switch1Idx] = volVars.moleFraction(FluidSystem::wPhaseIdx, FluidSystem::gCompIdx);
            }
            else if (volVars.saturation(FluidSystem::wPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // water phase disappears
                if (this->verbosity() > 1)
                    std::cout << "Water phase disappears at dof " << dofIdxGlobal
                              << ", coordinates: " << globalPos << ", sw: "
                              << volVars.saturation(FluidSystem::wPhaseIdx) << std::endl;
                newPhasePresence = Indices::gnPhaseOnly;

                priVars[Indices::switch1Idx] = volVars.moleFraction(FluidSystem::gPhaseIdx, FluidSystem::wCompIdx);
            }
            else if (volVars.saturation(FluidSystem::nPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // NAPL phase disappears
                if (this->verbosity() > 1)
                    std::cout << "NAPL phase disappears at dof " << dofIdxGlobal
                              << ", coordinates: " << globalPos << ", sn: "
                              << volVars.saturation(FluidSystem::nPhaseIdx) << std::endl;
                newPhasePresence = Indices::wgPhaseOnly;

                priVars[Indices::switch2Idx] = volVars.moleFraction(FluidSystem::gPhaseIdx, FluidSystem::nCompIdx);
            }
        }
        else if (phasePresence == Indices::wPhaseOnly)
        {
            bool gasPresent = false;
            bool nonwettingPresent = false;
            // calculate fractions of the partial pressures in the
            // hypothetical gas phase
            Scalar xwg = volVars.moleFraction(FluidSystem::gPhaseIdx, FluidSystem::wCompIdx);
            Scalar xgg = volVars.moleFraction(FluidSystem::gPhaseIdx, FluidSystem::gCompIdx);
            Scalar xng = volVars.moleFraction(FluidSystem::gPhaseIdx, FluidSystem::nCompIdx);
            /* take care:
               for xgg in case Indices::wPhaseOnly we compute xgg=henry_air*x2w
               for xwg in case Indices::wPhaseOnly we compute xwg=pwsat
               for xng in case Indices::wPhaseOnly we compute xng=henry_NAPL*x1w
            */

            Scalar xgMax = 1.0;
            if (xwg + xgg + xng > xgMax)
                wouldSwitch = true;
            if (this->wasSwitched_[dofIdxGlobal])
                xgMax *= 1.02;

            // if the sum of the mole fractions would be larger than
            // 100%, gas phase appears
            if (xwg + xgg + xng > xgMax)
            {
                // gas phase appears
                if (this->verbosity() > 1)
                    std::cout << "gas phase appears at dof " << dofIdxGlobal
                              << ", coordinates: " << globalPos << ", xwg + xgg + xng: "
                              << xwg + xgg + xng << std::endl;
                gasPresent = true;
            }

            // calculate fractions in the hypothetical NAPL phase
            Scalar xnn = volVars.moleFraction(FluidSystem::nPhaseIdx, FluidSystem::nCompIdx);
            /* take care:
               for xnn in case Indices::wPhaseOnly we compute xnn=henry_mesitylene*x1w,
               where a hypothetical gas pressure is assumed for the Henry
               x0n is set to NULL  (all NAPL phase is dirty)
               x2n is set to NULL  (all NAPL phase is dirty)
            */

            Scalar xnMax = 1.0;
            if (xnn > xnMax)
                wouldSwitch = true;
            if (this->wasSwitched_[dofIdxGlobal])
                xnMax *= 1.02;

            // if the sum of the hypothetical mole fractions would be larger than
            // 100%, NAPL phase appears
            if (xnn > xnMax)
            {
                // NAPL phase appears
                if (this->verbosity() > 1)
                    std::cout << "NAPL phase appears at dof " << dofIdxGlobal
                              << ", coordinates: " << globalPos << ", xnn: "
                              << xnn << std::endl;
                nonwettingPresent = true;
            }

            if (gasPresent && !nonwettingPresent)
            {
                newPhasePresence = Indices::wgPhaseOnly;
                priVars[Indices::switch1Idx] = 0.9999;
                priVars[Indices::switch2Idx] = 0.0001;
            }
            else if (gasPresent && nonwettingPresent)
            {
                newPhasePresence = Indices::threePhases;
                priVars[Indices::switch1Idx] = 0.9999;
                priVars[Indices::switch2Idx] = 0.0001;
            }
            else if (!gasPresent && nonwettingPresent)
            {
                newPhasePresence = Indices::wnPhaseOnly;
                priVars[Indices::switch1Idx] = volVars.moleFraction(FluidSystem::wPhaseIdx, FluidSystem::gCompIdx);
                priVars[Indices::switch2Idx] = 0.0001;
            }
        }
        else if (phasePresence == Indices::gnPhaseOnly)
        {
            bool nonwettingPresent = false;
            bool wettingPresent = false;

            Scalar Smin = 0.0;
            if (this->wasSwitched_[dofIdxGlobal])
                Smin = -0.01;

            if (volVars.saturation(FluidSystem::nPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // NAPL phase disappears
                if (this->verbosity() > 1)
                    std::cout << "NAPL phase disappears at dof " << dofIdxGlobal
                              << ", coordinates: " << globalPos << ", sn: "
                              << volVars.saturation(FluidSystem::nPhaseIdx) << std::endl;
                nonwettingPresent = true;
            }


            // calculate fractions of the hypothetical water phase
            Scalar xww = volVars.moleFraction(FluidSystem::wPhaseIdx, FluidSystem::wCompIdx);
            /*
              take care:, xww, if no water is present, then take xww=xwg*pg/pwsat .
              If this is larger than 1, then water appears
            */
            Scalar xwMax = 1.0;
            if (xww > xwMax)
                wouldSwitch = true;
            if (this->wasSwitched_[dofIdxGlobal])
                xwMax *= 1.02;

            // if the sum of the mole fractions would be larger than
            // 100%, gas phase appears
            if (xww > xwMax)
            {
                // water phase appears
                if (this->verbosity() > 1)
                    std::cout << "water phase appears at dof " << dofIdxGlobal
                              << ", coordinates: " << globalPos << ", xww=xwg*pg/pwsat : "
                              << xww << std::endl;
                wettingPresent = true;
            }

            if (wettingPresent && !nonwettingPresent)
            {
                newPhasePresence = Indices::threePhases;
                priVars[Indices::switch1Idx] = 0.0001;
                priVars[Indices::switch2Idx] = volVars.saturation(FluidSystem::nPhaseIdx);
            }
            else if (wettingPresent && nonwettingPresent)
            {
                newPhasePresence = Indices::wgPhaseOnly;
                priVars[Indices::switch1Idx] = 0.0001;
                priVars[Indices::switch2Idx] = volVars.moleFraction(FluidSystem::gPhaseIdx, FluidSystem::nCompIdx);
            }
            else if (!wettingPresent && nonwettingPresent)
            {
                newPhasePresence = Indices::gPhaseOnly;
                priVars[Indices::switch1Idx] = volVars.moleFraction(FluidSystem::gPhaseIdx, FluidSystem::wCompIdx);
                priVars[Indices::switch2Idx] = volVars.moleFraction(FluidSystem::gPhaseIdx, FluidSystem::nCompIdx);
            }
        }
        else if (phasePresence == Indices::wnPhaseOnly)
        {
            bool nonwettingPresent = false;
            bool gasPresent = false;

            Scalar Smin = 0.0;
            if (this->wasSwitched_[dofIdxGlobal])
                Smin = -0.01;

            if (volVars.saturation(FluidSystem::nPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // NAPL phase disappears
                if (this->verbosity() > 1)
                    std::cout << "NAPL phase disappears at dof " << dofIdxGlobal
                              << ", coordinates: " << globalPos << ", sn: "
                              << volVars.saturation(FluidSystem::nPhaseIdx) << std::endl;
                nonwettingPresent = true;
            }

            // calculate fractions of the partial pressures in the
            // hypothetical gas phase
            Scalar xwg = volVars.moleFraction(FluidSystem::gPhaseIdx, FluidSystem::wCompIdx);
            Scalar xgg = volVars.moleFraction(FluidSystem::gPhaseIdx, FluidSystem::gCompIdx);
            Scalar xng = volVars.moleFraction(FluidSystem::gPhaseIdx, FluidSystem::nCompIdx);
            /* take care:
               for xgg in case Indices::wPhaseOnly we compute xgg=henry_air*x2w
               for xwg in case Indices::wPhaseOnly we compute xwg=pwsat
               for xng in case Indices::wPhaseOnly we compute xng=henry_NAPL*x1w
            */
            Scalar xgMax = 1.0;
            if (xwg + xgg + xng > xgMax)
                wouldSwitch = true;
            if (this->wasSwitched_[dofIdxGlobal])
                xgMax *= 1.02;

            // if the sum of the mole fractions would be larger than
            // 100%, gas phase appears
            if (xwg + xgg + xng > xgMax)
            {
                // gas phase appears
                if (this->verbosity() > 1)
                    std::cout << "gas phase appears at dof " << dofIdxGlobal
                              << ", coordinates: " << globalPos << ", xwg + xgg + xng: "
                              << xwg + xgg + xng << std::endl;
                gasPresent = true;
            }

            if (gasPresent && !nonwettingPresent)
            {
                newPhasePresence = Indices::threePhases;
                priVars[Indices::switch1Idx] = volVars.saturation(FluidSystem::wPhaseIdx);
                priVars[Indices::switch2Idx] = volVars.saturation(FluidSystem::nPhaseIdx);
            }
            else if (gasPresent && nonwettingPresent)
            {
                newPhasePresence = Indices::wgPhaseOnly;
                priVars[Indices::switch1Idx] = volVars.saturation(FluidSystem::wPhaseIdx);
                priVars[Indices::switch2Idx] = volVars.moleFraction(FluidSystem::gPhaseIdx, FluidSystem::nCompIdx);
            }
            else if (!gasPresent && nonwettingPresent)
            {
                newPhasePresence = Indices::wPhaseOnly;
                priVars[Indices::switch1Idx] = volVars.moleFraction(FluidSystem::wPhaseIdx, FluidSystem::gCompIdx);
                priVars[Indices::switch2Idx] = volVars.moleFraction(FluidSystem::wPhaseIdx, FluidSystem::nCompIdx);
            }
        }
        else if (phasePresence == Indices::gPhaseOnly)
        {
            bool nonwettingPresent = false;
            bool wettingPresent = false;

            // calculate fractions in the hypothetical NAPL phase
            Scalar xnn = volVars.moleFraction(FluidSystem::nPhaseIdx, FluidSystem::nCompIdx);
            /*
              take care:, xnn, if no NAPL phase is there, take xnn=xng*pg/pcsat
              if this is larger than 1, then NAPL appears
            */

            Scalar xnMax = 1.0;
            if (xnn > xnMax)
                wouldSwitch = true;
            if (this->wasSwitched_[dofIdxGlobal])
                xnMax *= 1.02;

            // if the sum of the hypothetical mole fraction would be larger than
            // 100%, NAPL phase appears
            if (xnn > xnMax)
            {
                // NAPL phase appears
                if (this->verbosity() > 1)
                    std::cout << "NAPL phase appears at dof " << dofIdxGlobal
                              << ", coordinates: " << globalPos << ", xnn: "
                              << xnn << std::endl;
                nonwettingPresent = true;
            }
            // calculate fractions of the hypothetical water phase
            Scalar xww = volVars.moleFraction(FluidSystem::wPhaseIdx, FluidSystem::wCompIdx);
            /*
              take care:, xww, if no water is present, then take xww=xwg*pg/pwsat .
              If this is larger than 1, then water appears
            */
            Scalar xwMax = 1.0;
            if (xww > xwMax)
                wouldSwitch = true;
            if (this->wasSwitched_[dofIdxGlobal])
                xwMax *= 1.02;

            // if the sum of the mole fractions would be larger than
            // 100%, gas phase appears
            if (xww > xwMax)
            {
                // water phase appears
                if (this->verbosity() > 1)
                    std::cout << "water phase appears at dof " << dofIdxGlobal
                              << ", coordinates: " << globalPos << ", xww=xwg*pg/pwsat : "
                              << xww << std::endl;
                wettingPresent = true;
            }
            if (wettingPresent && !nonwettingPresent)
            {
                newPhasePresence = Indices::wgPhaseOnly;
                priVars[Indices::switch1Idx] = 0.0001;
                priVars[Indices::switch2Idx] = volVars.moleFraction(FluidSystem::gPhaseIdx, FluidSystem::nCompIdx);
            }
            else if (wettingPresent && nonwettingPresent)
            {
                newPhasePresence = Indices::threePhases;
                priVars[Indices::switch1Idx] = 0.0001;
                priVars[Indices::switch2Idx] = 0.0001;
            }
            else if (!wettingPresent && nonwettingPresent)
            {
                newPhasePresence = Indices::gnPhaseOnly;
                priVars[Indices::switch1Idx] = volVars.moleFraction(FluidSystem::gPhaseIdx, FluidSystem::wCompIdx);
                priVars[Indices::switch2Idx] = 0.0001;
            }
        }
        else if (phasePresence == Indices::wgPhaseOnly)
        {
            bool nonwettingPresent = false;
            bool gasPresent = false;
            bool wettingPresent = false;

            // get the fractions in the hypothetical NAPL phase
            Scalar xnn = volVars.moleFraction(FluidSystem::nPhaseIdx, FluidSystem::nCompIdx);

            // take care: if the NAPL phase is not present, take
            // xnn=xng*pg/pcsat if this is larger than 1, then NAPL
            // appears
            Scalar xnMax = 1.0;
            if (xnn > xnMax)
                wouldSwitch = true;
            if (this->wasSwitched_[dofIdxGlobal])
                xnMax *= 1.02;

            // if the sum of the hypothetical mole fraction would be larger than
            // 100%, NAPL phase appears
            if (xnn > xnMax)
            {
                // NAPL phase appears
                if (this->verbosity() > 1)
                    std::cout << "NAPL phase appears at dof " << dofIdxGlobal
                              << ", coordinates: " << globalPos << ", xnn: "
                              << xnn << std::endl;
                nonwettingPresent = true;
            }

            Scalar Smin = -1.e-6;
            if (this->wasSwitched_[dofIdxGlobal])
                Smin = -0.01;

            if (volVars.saturation(FluidSystem::gPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // gas phase disappears
                if (this->verbosity() > 1)
                    std::cout << "Gas phase disappears at dof " << dofIdxGlobal
                              << ", coordinates: " << globalPos << ", sg: "
                              << volVars.saturation(FluidSystem::gPhaseIdx) << std::endl;
                gasPresent = true;
            }

            Smin = 0.0;
            if (this->wasSwitched_[dofIdxGlobal])
                Smin = -0.01;

            if (volVars.saturation(FluidSystem::wPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // gas phase disappears
                if (this->verbosity() > 1)
                    std::cout << "Water phase disappears at dof " << dofIdxGlobal
                              << ", coordinates: " << globalPos << ", sw: "
                              << volVars.saturation(FluidSystem::wPhaseIdx) << std::endl;
                wettingPresent = true;
            }

            if (!gasPresent && nonwettingPresent && wettingPresent)
            {
                newPhasePresence = Indices::gnPhaseOnly;
                priVars[Indices::switch1Idx] = volVars.moleFraction(FluidSystem::gPhaseIdx, FluidSystem::wCompIdx);
                priVars[Indices::switch2Idx] = 0.0001;
            }
            else if (!gasPresent && nonwettingPresent && !wettingPresent)
            {
                newPhasePresence = Indices::threePhases;
                priVars[Indices::switch1Idx] = volVars.saturation(FluidSystem::wPhaseIdx);
                priVars[Indices::switch2Idx] = 0.0;
            }
            else if (gasPresent && !nonwettingPresent && !wettingPresent)
            {
                newPhasePresence = Indices::wPhaseOnly;
                priVars[Indices::switch1Idx] = volVars.moleFraction(FluidSystem::wPhaseIdx, FluidSystem::gCompIdx);
                priVars[Indices::switch2Idx] = volVars.moleFraction(FluidSystem::wPhaseIdx, FluidSystem::nCompIdx);
            }
            else if (!gasPresent && !nonwettingPresent && wettingPresent)
            {
                newPhasePresence = Indices::gPhaseOnly;
                priVars[Indices::switch1Idx] = volVars.moleFraction(FluidSystem::gPhaseIdx, FluidSystem::wCompIdx);
                priVars[Indices::switch2Idx] = volVars.moleFraction(FluidSystem::gPhaseIdx, FluidSystem::nCompIdx);
            }
        }

        priVars.setState(newPhasePresence);
        this->wasSwitched_[dofIdxGlobal] = wouldSwitch;
        return phasePresence != newPhasePresence;
    }
};

} // end namespace dumux

#endif
