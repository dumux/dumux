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
 * \ingroup ThreePWaterOilModel
 * \brief The primary variable switch for the 3p3c model
 */
#ifndef DUMUX_3P2CNI_PRIMARY_VARIABLE_SWITCH_HH
#define DUMUX_3P2CNI_PRIMARY_VARIABLE_SWITCH_HH

#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/compositional/primaryvariableswitch.hh>

namespace Dumux {

/*!
 * \ingroup ThreePWaterOilModel
 * \brief The primary variable switch controlling the phase presence state variable
 */
template<class TypeTag>
class ThreePWaterOilPrimaryVariableSwitch
: public PrimaryVariableSwitch<ThreePWaterOilPrimaryVariableSwitch<TypeTag>>
{
    using ParentType = PrimaryVariableSwitch<ThreePWaterOilPrimaryVariableSwitch<TypeTag>>;
    friend ParentType;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using GlobalPosition = Dune::FieldVector<Scalar, GridView::dimensionworld>;

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

    enum {
        switch1Idx = Indices::switch1Idx,
        switch2Idx = Indices::switch2Idx,
        pressureIdx = Indices::pressureIdx,

        wPhaseIdx = FluidSystem::wPhaseIdx,
        nPhaseIdx = FluidSystem::nPhaseIdx,
        gPhaseIdx = FluidSystem::gPhaseIdx,

        wCompIdx = FluidSystem::wCompIdx,
        nCompIdx = FluidSystem::nCompIdx,

        threePhases = Indices::threePhases,
        wPhaseOnly  = Indices::wPhaseOnly,
        gnPhaseOnly = Indices::gnPhaseOnly,
        wnPhaseOnly = Indices::wnPhaseOnly,
        gPhaseOnly  = Indices::gPhaseOnly,
        wgPhaseOnly = Indices::wgPhaseOnly
    };

public:
    using ParentType::ParentType;

protected:

    // perform variable switch at a degree of freedom location
    bool update_(PrimaryVariables& priVars,
                 const VolumeVariables& volVars,
                 IndexType dofIdxGlobal,
                 const GlobalPosition& globalPos)
    {
        // evaluate primary variable switch
        bool wouldSwitch = false;
        auto phasePresence = priVars.state();
        int newPhasePresence = phasePresence;

        bool onlyGasPhaseCanDisappear = GET_PROP_VALUE(TypeTag, OnlyGasPhaseCanDisappear);

        if(onlyGasPhaseCanDisappear)
        {
            // check if a primary var switch is necessary
            if (phasePresence == threePhases)
            {
                Scalar Smin = 0;
                if (this->wasSwitched_[dofIdxGlobal])
                    Smin = -0.01;

                if (volVars.saturation(gPhaseIdx) <= Smin)
                {
                    wouldSwitch = true;
                    // gas phase disappears
                    if (this->verbosity() > 1)
                        std::cout << "Gas phase disappears at dof " << dofIdxGlobal
                                << ", coordinates: " << globalPos << ", sg: "
                                << volVars.saturation(gPhaseIdx) << std::endl;
                    newPhasePresence = wnPhaseOnly;

                    priVars[pressureIdx] = volVars.fluidState().pressure(wPhaseIdx);
                    priVars[switch1Idx] = volVars.fluidState().temperature();
                }
            }
            else if (phasePresence == wnPhaseOnly)
            {
                bool gasFlag = 0;

                // calculate fractions of the partial pressures in the
                // hypothetical gas phase
                Scalar xwg = volVars.fluidState().moleFraction(gPhaseIdx, wCompIdx);
                Scalar xng = volVars.fluidState().moleFraction(gPhaseIdx, nCompIdx);

                Scalar xgMax = 1.0;
                if (xwg + xng > xgMax)
                    wouldSwitch = true;
                if (this->wasSwitched_[dofIdxGlobal])
                    xgMax *= 1.02;

                // if the sum of the mole fractions would be larger than
                // 100%, gas phase appears
                if (xwg + xng > xgMax)
                {
                    // gas phase appears
                    if (this->verbosity() > 1)
                        std::cout << "gas phase appears at dof " << dofIdxGlobal
                                << ", coordinates: " << globalPos << ", xwg + xng: "
                                << xwg + xng << std::endl;
                    gasFlag = 1;
                }

                if (gasFlag == 1)
                {
                    newPhasePresence = threePhases;
                    priVars[pressureIdx] = volVars.pressure(gPhaseIdx);
                    priVars[switch1Idx] = volVars.saturation(wPhaseIdx);
                }
            }
        }

        else
        {
            // check if a primary var switch is necessary
            if (phasePresence == threePhases)
            {
                Scalar Smin = 0;
                if (this->wasSwitched_[dofIdxGlobal])
                    Smin = -0.01;

                if (volVars.saturation(gPhaseIdx) <= Smin)
                {
                    wouldSwitch = true;
                    // gas phase disappears
                    if (this->verbosity() > 1)
                        std::cout << "Gas phase disappears at dof " << dofIdxGlobal
                                << ", coordinates: " << globalPos << ", sg: "
                                << volVars.saturation(gPhaseIdx) << std::endl;
                    newPhasePresence = wnPhaseOnly;

                    priVars[pressureIdx] = volVars.fluidState().pressure(wPhaseIdx);
                    priVars[switch1Idx] = volVars.fluidState().temperature();
                }
                else if (volVars.saturation(wPhaseIdx) <= Smin)
                {
                    wouldSwitch = true;
                    // water phase disappears
                    if (this->verbosity() > 1)
                        std::cout << "Water phase disappears at dof " << dofIdxGlobal
                                << ", coordinates: " << globalPos << ", sw: "
                                << volVars.saturation(wPhaseIdx) << std::endl;
                    newPhasePresence = gnPhaseOnly;

                    priVars[switch1Idx] = volVars.fluidState().saturation(nPhaseIdx);
                    priVars[switch2Idx] = volVars.fluidState().moleFraction(nPhaseIdx, wCompIdx);
                }
                else if (volVars.saturation(nPhaseIdx) <= Smin)
                {
                    wouldSwitch = true;
                    // NAPL phase disappears
                    if (this->verbosity() > 1)
                        std::cout << "NAPL phase disappears at dof " << dofIdxGlobal
                                << ", coordinates: " << globalPos << ", sn: "
                                << volVars.saturation(nPhaseIdx) << std::endl;
                    newPhasePresence = wgPhaseOnly;

                    priVars[switch2Idx] = volVars.fluidState().moleFraction(wPhaseIdx, nCompIdx);
                }
            }
            else if (phasePresence == wPhaseOnly)
            {
                bool gasFlag = 0;
                bool nonwettingFlag = 0;
                // calculate fractions of the partial pressures in the
                // hypothetical gas phase
                Scalar xwg = volVars.fluidState().moleFraction(gPhaseIdx, wCompIdx);
                Scalar xng = volVars.fluidState().moleFraction(gPhaseIdx, nCompIdx);

                Scalar xgMax = 1.0;
                if (xwg + xng > xgMax)
                    wouldSwitch = true;
                if (this->wasSwitched_[dofIdxGlobal])
                    xgMax *= 1.02;

                // if the sum of the mole fractions would be larger than
                // 100%, gas phase appears
                if (xwg + xng > xgMax)
                {
                    // gas phase appears
                    if (this->verbosity() > 1)
                        std::cout << "gas phase appears at dof " << dofIdxGlobal
                                << ", coordinates: " << globalPos << ", xwg + xng: "
                                << xwg + xng << std::endl;
                    gasFlag = 1;
                }

                // calculate fractions in the hypothetical NAPL phase
                Scalar xnn = volVars.fluidState().moleFraction(nPhaseIdx, nCompIdx);
                /* take care:
                for xnn in case wPhaseOnly we compute xnn=henry_mesitylene*x1w,
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
                    nonwettingFlag = 1;
                }

                if ((gasFlag == 1) && (nonwettingFlag == 0))
                {
                    newPhasePresence = wgPhaseOnly;
                    priVars[switch1Idx] = 0.9999;
                    priVars[pressureIdx] = volVars.fluidState().pressure(gPhaseIdx);
                }
                else if ((gasFlag == 1) && (nonwettingFlag == 1))
                {
                    newPhasePresence = threePhases;
                    priVars[pressureIdx] = volVars.fluidState().pressure(gPhaseIdx);
                    priVars[switch1Idx] = 0.999;
                }
                else if ((gasFlag == 0) && (nonwettingFlag == 1))
                {
                    newPhasePresence = wnPhaseOnly;
                    priVars[switch2Idx] = 0.0001;
                }
            }
            else if (phasePresence == gnPhaseOnly)
            {
                bool nonwettingFlag = 0;
                bool wettingFlag = 0;

                Scalar Smin = 0.0;
                if (this->wasSwitched_[dofIdxGlobal])
                    Smin = -0.01;

                if (volVars.saturation(nPhaseIdx) <= Smin)
                {
                    wouldSwitch = true;
                    // NAPL phase disappears
                    if (this->verbosity() > 1)
                        std::cout << "NAPL phase disappears at dof " << dofIdxGlobal
                                << ", coordinates: " << globalPos << ", sn: "
                                << volVars.saturation(nPhaseIdx) << std::endl;
                    nonwettingFlag = 1;
                }


                // calculate fractions of the hypothetical water phase
                Scalar xww = volVars.fluidState().moleFraction(wPhaseIdx, wCompIdx);
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
                // 100%, water phase appears
                if (xww > xwMax)
                {
                    // water phase appears
                    if (this->verbosity() > 1)
                        std::cout << "water phase appears at dof " << dofIdxGlobal
                                << ", coordinates: " << globalPos << ", xww=xwg*pg/pwsat : "
                                << xww << std::endl;
                    wettingFlag = 1;
                }

                if ((wettingFlag == 1) && (nonwettingFlag == 0))
                {
                    newPhasePresence = threePhases;
                    priVars[switch1Idx] = 0.0001;
                    priVars[switch2Idx] = volVars.saturation(nPhaseIdx);
                }
                else if ((wettingFlag == 1) && (nonwettingFlag == 1))
                {
                    newPhasePresence = wgPhaseOnly;
                    priVars[switch1Idx] = 0.0001;
                    priVars[switch2Idx] = volVars.fluidState().moleFraction(wPhaseIdx, nCompIdx);
                }
                else if ((wettingFlag == 0) && (nonwettingFlag == 1))
                {
                    newPhasePresence = gPhaseOnly;
                    priVars[switch1Idx] = volVars.fluidState().temperature();
                    priVars[switch2Idx] = volVars.fluidState().moleFraction(gPhaseIdx, nCompIdx);
                }
            }
            else if (phasePresence == wnPhaseOnly)
            {
                bool nonwettingFlag = 0;
                bool gasFlag = 0;

                Scalar Smin = 0.0;
                if (this->wasSwitched_[dofIdxGlobal])
                    Smin = -0.01;

                if (volVars.saturation(nPhaseIdx) <= Smin)
                {
                    wouldSwitch = true;
                    // NAPL phase disappears
                    if (this->verbosity() > 1)
                        std::cout << "NAPL phase disappears at dof " << dofIdxGlobal
                                << ", coordinates: " << globalPos << ", sn: "
                                << volVars.saturation(nPhaseIdx) << std::endl;
                    nonwettingFlag = 1;
                }

                // calculate fractions of the partial pressures in the
                // hypothetical gas phase
                Scalar xwg = volVars.fluidState().moleFraction(gPhaseIdx, wCompIdx);
                Scalar xng = volVars.fluidState().moleFraction(gPhaseIdx, nCompIdx);

                Scalar xgMax = 1.0;
                if (xwg + xng > xgMax)
                    wouldSwitch = true;
                if (this->wasSwitched_[dofIdxGlobal])
                    xgMax *= 1.02;

                // if the sum of the mole fractions would be larger than
                // 100%, gas phase appears
                if (xwg + xng > xgMax)
                {
                    // gas phase appears
                    if (this->verbosity() > 1)
                        std::cout << "gas phase appears at dof " << dofIdxGlobal
                                << ", coordinates: " << globalPos << ", xwg + xng: "
                                << xwg + xng << std::endl;
                    gasFlag = 1;
                }

                if ((gasFlag == 1) && (nonwettingFlag == 0))
                {
                    newPhasePresence = threePhases;
                    priVars[pressureIdx] = volVars.pressure(gPhaseIdx);
                    priVars[switch1Idx] = volVars.saturation(wPhaseIdx);
                }
                else if ((gasFlag == 1) && (nonwettingFlag == 1))
                {
                    newPhasePresence = wgPhaseOnly;
                    priVars[pressureIdx] = volVars.pressure(gPhaseIdx);
                    priVars[switch1Idx] = volVars.temperature();
                    priVars[switch2Idx] = volVars.fluidState().moleFraction(wPhaseIdx, nCompIdx);
                }
                else if ((gasFlag == 0) && (nonwettingFlag == 1))
                {
                    newPhasePresence = wPhaseOnly;
                    priVars[switch2Idx] = volVars.fluidState().moleFraction(wPhaseIdx, nCompIdx);
                }
            }
            else if (phasePresence == gPhaseOnly)
            {
                bool nonwettingFlag = 0;
                bool wettingFlag = 0;

                // calculate fractions in the hypothetical NAPL phase
                Scalar xnn = volVars.fluidState().moleFraction(nPhaseIdx, nCompIdx);
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
                    nonwettingFlag = 1;
                }
                // calculate fractions of the hypothetical water phase
                Scalar xww = volVars.fluidState().moleFraction(wPhaseIdx, wCompIdx);
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
                    wettingFlag = 1;
                }
                if ((wettingFlag == 1) && (nonwettingFlag == 0))
                {
                    newPhasePresence = wgPhaseOnly;
                    priVars[switch1Idx] = 0.0001;
                    priVars[switch2Idx] = volVars.fluidState().moleFraction(wPhaseIdx, nCompIdx);
                }
                else if ((wettingFlag == 1) && (nonwettingFlag == 1))
                {
                    newPhasePresence = threePhases;
                    priVars[switch1Idx] = 0.0001;
                    priVars[switch2Idx] = 0.0001;
                }
                else if ((wettingFlag == 0) && (nonwettingFlag == 1))
                {
                    newPhasePresence = gnPhaseOnly;
                    priVars[switch2Idx]
                        = volVars.fluidState().moleFraction(wPhaseIdx, nCompIdx);
                }
            }
            else if (phasePresence == wgPhaseOnly)
            {
                bool nonwettingFlag = 0;
                bool gasFlag = 0;
                bool wettingFlag = 0;

                // get the fractions in the hypothetical NAPL phase
                Scalar xnn = volVars.fluidState().moleFraction(nPhaseIdx, nCompIdx);

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
                    nonwettingFlag = 1;
                }

                Scalar Smin = -1.e-6;
                if (this->wasSwitched_[dofIdxGlobal])
                    Smin = -0.01;

                if (volVars.saturation(gPhaseIdx) <= Smin)
                {
                    wouldSwitch = true;
                    // gas phase disappears
                    if (this->verbosity() > 1)
                        std::cout << "Gas phase disappears at dof " << dofIdxGlobal
                                << ", coordinates: " << globalPos << ", sg: "
                                << volVars.saturation(gPhaseIdx) << std::endl;
                    gasFlag = 1;
                }

                Smin = 0.0;
                if (this->wasSwitched_[dofIdxGlobal])
                    Smin = -0.01;

                if (volVars.saturation(wPhaseIdx) <= Smin)
                {
                    wouldSwitch = true;
                    // gas phase disappears
                    if (this->verbosity() > 1)
                        std::cout << "Water phase disappears at dof " << dofIdxGlobal
                                << ", coordinates: " << globalPos << ", sw: "
                                << volVars.saturation(wPhaseIdx) << std::endl;
                    wettingFlag = 1;
                }

                if ((gasFlag == 0) && (nonwettingFlag == 1) && (wettingFlag == 1))
                {
                    newPhasePresence = gnPhaseOnly;
                    priVars[switch1Idx] = 0.0001;
                    priVars[switch2Idx] = volVars.fluidState().moleFraction(nPhaseIdx, wCompIdx);
    ;
                }
                else if ((gasFlag == 0) && (nonwettingFlag == 1) && (wettingFlag == 0))
                {
                    newPhasePresence = threePhases;
                    priVars[switch2Idx] = 0.0001;
                }
                else if ((gasFlag == 1) && (nonwettingFlag == 0) && (wettingFlag == 0))
                {
                    newPhasePresence = wPhaseOnly;
                    priVars[pressureIdx] = volVars.fluidState().pressure(wPhaseIdx);
                    priVars[switch1Idx] = volVars.fluidState().temperature();
                    priVars[switch2Idx] = volVars.fluidState().moleFraction(wPhaseIdx, nCompIdx);
                }
                else if ((gasFlag == 0) && (nonwettingFlag == 0) && (wettingFlag == 1))
                {
                    newPhasePresence = gPhaseOnly;
                    priVars[switch1Idx]
                        = volVars.fluidState().temperature();
                    priVars[switch2Idx]
                        = volVars.fluidState().moleFraction(gPhaseIdx, nCompIdx);
                }
            }
        }


        priVars.setState(newPhasePresence);
        this->wasSwitched_[dofIdxGlobal] = wouldSwitch;
        return phasePresence != newPhasePresence;
    }

};

} // end namespace dumux

#endif
