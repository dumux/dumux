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
 *
 * \brief The primary variable switch for the 2pnc model
 */
#ifndef DUMUX_2PNC_PRIMARY_VARIABLE_SWITCH_HH
#define DUMUX_2PNC_PRIMARY_VARIABLE_SWITCH_HH

#include <dumux/porousmediumflow/compositional/primaryvariableswitch.hh>

namespace Dumux
{
/*!
 * \ingroup TwoPNCModel
 * \brief The primary variable switch controlling the phase presence state variable
 */
template<class TypeTag>
class TwoPNCPrimaryVariableSwitch : public Dumux::PrimaryVariableSwitch<TypeTag>
{
    friend typename Dumux::PrimaryVariableSwitch<TypeTag>;
    using ParentType = Dumux::PrimaryVariableSwitch<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using GlobalPosition = Dune::FieldVector<Scalar, GridView::dimensionworld>;

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    static const int numComponents = GET_PROP_VALUE(TypeTag, NumComponents);
    static const int numMajorComponents = GET_PROP_VALUE(TypeTag, NumMajorComponents);

    enum {
        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx,

        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,

        wCompIdx = FluidSystem::wCompIdx,
        nCompIdx = FluidSystem::nCompIdx,

        wPhaseOnly = Indices::wPhaseOnly,
        nPhaseOnly = Indices::nPhaseOnly,
        bothPhases = Indices::bothPhases
    };

    enum {
            plSg = TwoPNCFormulation::plSg,
            pgSl = TwoPNCFormulation::pgSl,
            formulation = GET_PROP_VALUE(TypeTag, Formulation)
    };

protected:

    // perform variable switch at a degree of freedom location
    bool update_(PrimaryVariables& priVars,
                 const VolumeVariables& volVars,
                 IndexType dofIdxGlobal,
                 const GlobalPosition& globalPos)
    {
        // evaluate primary variable switch
        bool wouldSwitch = false;
        int phasePresence = this->phasePresence_[dofIdxGlobal];
        int newPhasePresence = phasePresence;

        //check if a primary variable switch is necessary
        if (phasePresence == bothPhases)
        {
            Scalar Smin = 0; //saturation threshold
            if (this->wasSwitched_[dofIdxGlobal])
                Smin = -0.01;

            //if saturation of liquid phase is smaller 0 switch
            if (volVars.saturation(wPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                //liquid phase has to disappear
                std::cout << "Liquid Phase disappears at vertex " << dofIdxGlobal
                            << ", coordinated: " << globalPos << ", Sl: "
                            << volVars.saturation(wPhaseIdx) << std::endl;
                newPhasePresence = nPhaseOnly;

                //switch not depending on formulation
                //switch "Sl" to "xgH20"
                priVars[switchIdx]
                        = volVars.moleFraction(nPhaseIdx, wCompIdx /*H2O*/);

                //switch all secondary components to mole fraction in gas phase
                for (int compIdx=numMajorComponents; compIdx<numComponents; ++compIdx)
                    priVars[compIdx] = volVars.moleFraction(nPhaseIdx,compIdx);
            }
            //if saturation of gas phase is smaller than 0 switch
            else if (volVars.saturation(nPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                //gas phase has to disappear
                std::cout << "Gas Phase disappears at vertex " << dofIdxGlobal
                            << ", coordinated: " << globalPos << ", Sg: "
                            << volVars.saturation(nPhaseIdx) << std::endl;
                newPhasePresence = wPhaseOnly;

                //switch "Sl" to "xlN2"
                priVars[switchIdx] = volVars.moleFraction(wPhaseIdx, nCompIdx /*N2*/);
            }
        }
        else if (phasePresence == nPhaseOnly)
        {
            Scalar xlmax = 1;
            Scalar sumxl = 0;
            //Calculate sum of mole fractions in the hypothetical liquid phase
            for (int compIdx = 0; compIdx < numComponents; compIdx++)
            {
                sumxl += volVars.moleFraction(wPhaseIdx, compIdx);
            }
            if (sumxl > xlmax)
                wouldSwitch = true;
            if (this->wasSwitched_[dofIdxGlobal])
                xlmax *=1.02;
            //liquid phase appears if sum is larger than one
            if (sumxl/*sum of mole fractions*/ > xlmax/*1*/)
            {
                std::cout << "Liquid Phase appears at vertex " << dofIdxGlobal
                        << ", coordinated: " << globalPos << ", sumxl: "
                        << sumxl << std::endl;
                newPhasePresence = bothPhases;

                //saturation of the liquid phase set to 0.0001 (if formulation pgSl and vice versa)
                if (formulation == pgSl)
                    priVars[switchIdx] = 0.0001;
                else if (formulation == plSg)
                    priVars[switchIdx] = 0.9999;

                //switch all secondary components back to liquid mole fraction
                for (int compIdx=numMajorComponents; compIdx<numComponents; ++compIdx)
                    priVars[compIdx] = volVars.moleFraction(wPhaseIdx,compIdx);
            }
        }
        else if (phasePresence == wPhaseOnly)
        {
            Scalar xgmax = 1;
            Scalar sumxg = 0;
            //Calculate sum of mole fractions in the hypothetical liquid phase
            for (int compIdx = 0; compIdx < numComponents; compIdx++)
            {
                sumxg += volVars.moleFraction(nPhaseIdx, compIdx);
            }
            if (sumxg > xgmax)
                wouldSwitch = true;
            if (this->wasSwitched_[dofIdxGlobal])
                xgmax *=1.02;
            //liquid phase appears if sum is larger than one
            if (sumxg > xgmax)
            {
                std::cout << "Gas Phase appears at vertex " << dofIdxGlobal
                        << ", coordinated: " << globalPos << ", sumxg: "
                        << sumxg << std::endl;
                newPhasePresence = bothPhases;
                //saturation of the liquid phase set to 0.9999 (if formulation pgSl and vice versa)
                if (formulation == pgSl)
                    priVars[switchIdx] = 0.9999;
                else if (formulation == plSg)
                    priVars[switchIdx] = 0.0001;

            }
        }


        this->phasePresence_[dofIdxGlobal] = newPhasePresence;
        this->wasSwitched_[dofIdxGlobal] = wouldSwitch;
        return phasePresence != newPhasePresence;
    }
};

} // end namespace Dumux

#endif
