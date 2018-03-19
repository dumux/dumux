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
 * \ingroup ThreePThreeCModel
 * \brief The primary variable switch for the 3p3c model
 */
#ifndef DUMUX_3P3C_PRIMARY_VARIABLE_SWITCH_HH
#define DUMUX_3P3C_PRIMARY_VARIABLE_SWITCH_HH

#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/compositional/primaryvariableswitch.hh>

namespace Dumux {

/*!
 * \ingroup ThreePThreeCModel
 * \brief The primary variable switch controlling the phase presence state variable
 */
template<class TypeTag>
class ThreePThreeCPrimaryVariableSwitch
: public PrimaryVariableSwitch<typename GET_PROP_TYPE(TypeTag, FVGridGeometry), ThreePThreeCPrimaryVariableSwitch<TypeTag>>
{
    using ParentType = PrimaryVariableSwitch<typename GET_PROP_TYPE(TypeTag, FVGridGeometry), ThreePThreeCPrimaryVariableSwitch<TypeTag>>;
    friend ParentType;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using GlobalPosition = Dune::FieldVector<Scalar, GridView::dimensionworld>;

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    enum {
        switch1Idx = Indices::switch1Idx,
        switch2Idx = Indices::switch2Idx,

        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        gPhaseIdx = Indices::gPhaseIdx,

        wCompIdx = Indices::wCompIdx,
        nCompIdx = Indices::nCompIdx,
        gCompIdx = Indices::gCompIdx,

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
        // evaluate if the primary variable switch would switch
        bool wouldSwitch = false;
        auto phasePresence = priVars.state();
        auto newPhasePresence = phasePresence;

        // check if a primary variable switch is necessary
        if (phasePresence == threePhases)
        {
            Scalar Smin = 0;
            if (this->wasSwitched_[dofIdxGlobal])
                Smin = -0.01;

            if (volVars.saturation(gPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // gas phase disappears
                std::cout << "Gas phase disappears at dof " << dofIdxGlobal
                          << ", coordinates: " << globalPos << ", sg: "
                          << volVars.saturation(gPhaseIdx) << std::endl;
                newPhasePresence = wnPhaseOnly;

                priVars[switch1Idx] = volVars.moleFraction(wPhaseIdx, gCompIdx);
            }
            else if (volVars.saturation(wPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // water phase disappears
                std::cout << "Water phase disappears at dof " << dofIdxGlobal
                          << ", coordinates: " << globalPos << ", sw: "
                          << volVars.saturation(wPhaseIdx) << std::endl;
                newPhasePresence = gnPhaseOnly;

                priVars[switch1Idx] = volVars.moleFraction(gPhaseIdx, wCompIdx);
            }
            else if (volVars.saturation(nPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // NAPL phase disappears
                std::cout << "NAPL phase disappears at dof " << dofIdxGlobal
                          << ", coordinates: " << globalPos << ", sn: "
                          << volVars.saturation(nPhaseIdx) << std::endl;
                newPhasePresence = wgPhaseOnly;

                priVars[switch2Idx] = volVars.moleFraction(gPhaseIdx, nCompIdx);
            }
        }
        else if (phasePresence == wPhaseOnly)
        {
            bool gasPresent = false;
            bool nonwettingPresent = false;
            // calculate fractions of the partial pressures in the
            // hypothetical gas phase
            Scalar xwg = volVars.moleFraction(gPhaseIdx, wCompIdx);
            Scalar xgg = volVars.moleFraction(gPhaseIdx, gCompIdx);
            Scalar xng = volVars.moleFraction(gPhaseIdx, nCompIdx);
            /* take care:
               for xgg in case wPhaseOnly we compute xgg=henry_air*x2w
               for xwg in case wPhaseOnly we compute xwg=pwsat
               for xng in case wPhaseOnly we compute xng=henry_NAPL*x1w
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
                std::cout << "gas phase appears at dof " << dofIdxGlobal
                          << ", coordinates: " << globalPos << ", xwg + xgg + xng: "
                          << xwg + xgg + xng << std::endl;
                gasPresent = true;
            }

            // calculate fractions in the hypothetical NAPL phase
            Scalar xnn = volVars.moleFraction(nPhaseIdx, nCompIdx);
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
                std::cout << "NAPL phase appears at dof " << dofIdxGlobal
                          << ", coordinates: " << globalPos << ", xnn: "
                          << xnn << std::endl;
                nonwettingPresent = true;
            }

            if (gasPresent && !nonwettingPresent)
            {
                newPhasePresence = wgPhaseOnly;
                priVars[switch1Idx] = 0.9999;
                priVars[switch2Idx] = 0.0001;
            }
            else if (gasPresent && nonwettingPresent)
            {
                newPhasePresence = threePhases;
                priVars[switch1Idx] = 0.9999;
                priVars[switch2Idx] = 0.0001;
            }
            else if (!gasPresent && nonwettingPresent)
            {
                newPhasePresence = wnPhaseOnly;
                priVars[switch1Idx] = volVars.moleFraction(wPhaseIdx, gCompIdx);
                priVars[switch2Idx] = 0.0001;
            }
        }
        else if (phasePresence == gnPhaseOnly)
        {
            bool nonwettingPresent = false;
            bool wettingPresent = false;

            Scalar Smin = 0.0;
            if (this->wasSwitched_[dofIdxGlobal])
                Smin = -0.01;

            if (volVars.saturation(nPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // NAPL phase disappears
                std::cout << "NAPL phase disappears at dof " << dofIdxGlobal
                          << ", coordinates: " << globalPos << ", sn: "
                          << volVars.saturation(nPhaseIdx) << std::endl;
                nonwettingPresent = true;
            }


            // calculate fractions of the hypothetical water phase
            Scalar xww = volVars.moleFraction(wPhaseIdx, wCompIdx);
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
                std::cout << "water phase appears at dof " << dofIdxGlobal
                          << ", coordinates: " << globalPos << ", xww=xwg*pg/pwsat : "
                          << xww << std::endl;
                wettingPresent = true;
            }

            if (wettingPresent && !nonwettingPresent)
            {
                newPhasePresence = threePhases;
                priVars[switch1Idx] = 0.0001;
                priVars[switch2Idx] = volVars.saturation(nPhaseIdx);
            }
            else if (wettingPresent && nonwettingPresent)
            {
                newPhasePresence = wgPhaseOnly;
                priVars[switch1Idx] = 0.0001;
                priVars[switch2Idx] = volVars.moleFraction(gPhaseIdx, nCompIdx);
            }
            else if (!wettingPresent && nonwettingPresent)
            {
                newPhasePresence = gPhaseOnly;
                priVars[switch1Idx] = volVars.moleFraction(gPhaseIdx, wCompIdx);
                priVars[switch2Idx] = volVars.moleFraction(gPhaseIdx, nCompIdx);
            }
        }
        else if (phasePresence == wnPhaseOnly)
        {
            bool nonwettingPresent = false;
            bool gasPresent = false;

            Scalar Smin = 0.0;
            if (this->wasSwitched_[dofIdxGlobal])
                Smin = -0.01;

            if (volVars.saturation(nPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // NAPL phase disappears
                std::cout << "NAPL phase disappears at dof " << dofIdxGlobal
                          << ", coordinates: " << globalPos << ", sn: "
                          << volVars.saturation(nPhaseIdx) << std::endl;
                nonwettingPresent = true;
            }

            // calculate fractions of the partial pressures in the
            // hypothetical gas phase
            Scalar xwg = volVars.moleFraction(gPhaseIdx, wCompIdx);
            Scalar xgg = volVars.moleFraction(gPhaseIdx, gCompIdx);
            Scalar xng = volVars.moleFraction(gPhaseIdx, nCompIdx);
            /* take care:
               for xgg in case wPhaseOnly we compute xgg=henry_air*x2w
               for xwg in case wPhaseOnly we compute xwg=pwsat
               for xng in case wPhaseOnly we compute xng=henry_NAPL*x1w
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
                std::cout << "gas phase appears at dof " << dofIdxGlobal
                          << ", coordinates: " << globalPos << ", xwg + xgg + xng: "
                          << xwg + xgg + xng << std::endl;
                gasPresent = true;
            }

            if (gasPresent && !nonwettingPresent)
            {
                newPhasePresence = threePhases;
                priVars[switch1Idx] = volVars.saturation(wPhaseIdx);
                priVars[switch2Idx] = volVars.saturation(nPhaseIdx);
            }
            else if (gasPresent && nonwettingPresent)
            {
                newPhasePresence = wgPhaseOnly;
                priVars[switch1Idx] = volVars.saturation(wPhaseIdx);
                priVars[switch2Idx] = volVars.moleFraction(gPhaseIdx, nCompIdx);
            }
            else if (!gasPresent && nonwettingPresent)
            {
                newPhasePresence = wPhaseOnly;
                priVars[switch1Idx] = volVars.moleFraction(wPhaseIdx, gCompIdx);
                priVars[switch2Idx] = volVars.moleFraction(wPhaseIdx, nCompIdx);
            }
        }
        else if (phasePresence == gPhaseOnly)
        {
            bool nonwettingPresent = false;
            bool wettingPresent = false;

            // calculate fractions in the hypothetical NAPL phase
            Scalar xnn = volVars.moleFraction(nPhaseIdx, nCompIdx);
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
                std::cout << "NAPL phase appears at dof " << dofIdxGlobal
                          << ", coordinates: " << globalPos << ", xnn: "
                          << xnn << std::endl;
                nonwettingPresent = true;
            }
            // calculate fractions of the hypothetical water phase
            Scalar xww = volVars.moleFraction(wPhaseIdx, wCompIdx);
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
                std::cout << "water phase appears at dof " << dofIdxGlobal
                          << ", coordinates: " << globalPos << ", xww=xwg*pg/pwsat : "
                          << xww << std::endl;
                wettingPresent = true;
            }
            if (wettingPresent && !nonwettingPresent)
            {
                newPhasePresence = wgPhaseOnly;
                priVars[switch1Idx] = 0.0001;
                priVars[switch2Idx] = volVars.moleFraction(gPhaseIdx, nCompIdx);
            }
            else if (wettingPresent && nonwettingPresent)
            {
                newPhasePresence = threePhases;
                priVars[switch1Idx] = 0.0001;
                priVars[switch2Idx] = 0.0001;
            }
            else if (!wettingPresent && nonwettingPresent)
            {
                newPhasePresence = gnPhaseOnly;
                priVars[switch1Idx] = volVars.moleFraction(gPhaseIdx, wCompIdx);
                priVars[switch2Idx] = 0.0001;
            }
        }
        else if (phasePresence == wgPhaseOnly)
        {
            bool nonwettingPresent = false;
            bool gasPresent = false;
            bool wettingPresent = false;

            // get the fractions in the hypothetical NAPL phase
            Scalar xnn = volVars.moleFraction(nPhaseIdx, nCompIdx);

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
                std::cout << "NAPL phase appears at dof " << dofIdxGlobal
                          << ", coordinates: " << globalPos << ", xnn: "
                          << xnn << std::endl;
                nonwettingPresent = true;
            }

            Scalar Smin = -1.e-6;
            if (this->wasSwitched_[dofIdxGlobal])
                Smin = -0.01;

            if (volVars.saturation(gPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // gas phase disappears
                std::cout << "Gas phase disappears at dof " << dofIdxGlobal
                          << ", coordinates: " << globalPos << ", sg: "
                          << volVars.saturation(gPhaseIdx) << std::endl;
                gasPresent = true;
            }

            Smin = 0.0;
            if (this->wasSwitched_[dofIdxGlobal])
                Smin = -0.01;

            if (volVars.saturation(wPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // gas phase disappears
                std::cout << "Water phase disappears at dof " << dofIdxGlobal
                          << ", coordinates: " << globalPos << ", sw: "
                          << volVars.saturation(wPhaseIdx) << std::endl;
                wettingPresent = true;
            }

            if (!gasPresent && nonwettingPresent && wettingPresent)
            {
                newPhasePresence = gnPhaseOnly;
                priVars[switch1Idx] = volVars.moleFraction(gPhaseIdx, wCompIdx);
                priVars[switch2Idx] = 0.0001;
            }
            else if (!gasPresent && nonwettingPresent && !wettingPresent)
            {
                newPhasePresence = threePhases;
                priVars[switch1Idx] = volVars.saturation(wPhaseIdx);
                priVars[switch2Idx] = 0.0;
            }
            else if (gasPresent && !nonwettingPresent && !wettingPresent)
            {
                newPhasePresence = wPhaseOnly;
                priVars[switch1Idx] = volVars.moleFraction(wPhaseIdx, gCompIdx);
                priVars[switch2Idx] = volVars.moleFraction(wPhaseIdx, nCompIdx);
            }
            else if (!gasPresent && !nonwettingPresent && wettingPresent)
            {
                newPhasePresence = gPhaseOnly;
                priVars[switch1Idx] = volVars.moleFraction(gPhaseIdx, wCompIdx);
                priVars[switch2Idx] = volVars.moleFraction(gPhaseIdx, nCompIdx);
            }
        }

        priVars.setState(newPhasePresence);
        this->wasSwitched_[dofIdxGlobal] = wouldSwitch;
        return phasePresence != newPhasePresence;
    }
};

} // end namespace dumux

#endif
