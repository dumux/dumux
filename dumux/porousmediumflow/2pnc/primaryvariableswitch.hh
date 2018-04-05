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
 * \ingroup TwoPNCModel
 * \brief The primary variable switch for the 2pnc model
 */
#ifndef DUMUX_2PNC_PRIMARY_VARIABLE_SWITCH_HH
#define DUMUX_2PNC_PRIMARY_VARIABLE_SWITCH_HH

#include <dumux/porousmediumflow/compositional/primaryvariableswitch.hh>
#include <dumux/porousmediumflow/2p/formulation.hh>

namespace Dumux {

/*!
 * \ingroup TwoPNCModel
 * \brief The primary variable switch controlling the phase presence state variable
 */
template<class TypeTag>
class TwoPNCPrimaryVariableSwitch
: public PrimaryVariableSwitch<TwoPNCPrimaryVariableSwitch<TypeTag>>
{
    using ParentType = PrimaryVariableSwitch<TwoPNCPrimaryVariableSwitch<TypeTag>>;
    friend ParentType;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using GlobalPosition = Dune::FieldVector<Scalar, GridView::dimensionworld>;

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;

    static const int numComponents = GET_PROP_TYPE(TypeTag, ModelTraits)::numComponents();
    static const int numMajorComponents = GET_PROP_TYPE(TypeTag, ModelTraits)::numMajorComponents();

    enum {

        switchIdx = Indices::switchIdx,

        wPhaseIdx = FluidSystem::wPhaseIdx,
        nPhaseIdx = FluidSystem::nPhaseIdx,

        wCompIdx = FluidSystem::wCompIdx,
        nCompIdx = FluidSystem::nCompIdx,

        wPhaseOnly = Indices::wPhaseOnly,
        nPhaseOnly = Indices::nPhaseOnly,
        bothPhases = Indices::bothPhases
    };

    static constexpr auto formulation = GET_PROP_TYPE(TypeTag, ModelTraits)::priVarFormulation();
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
        int phasePresence = priVars.state();
        int newPhasePresence = phasePresence;

        //check if a primary variable switch is necessary
        if (phasePresence == bothPhases)
        {
            Scalar Smin = 0; //saturation threshold
            if (this->wasSwitched_[dofIdxGlobal])
                Smin = -0.01;

            //if saturation of wetting phase is smaller 0 switch
            if (volVars.saturation(wPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                //wetting phase has to disappear
                std::cout << "Wetting Phase disappears at vertex " << dofIdxGlobal
                            << ", coordinated: " << globalPos << ", Sw: "
                            << volVars.saturation(wPhaseIdx) << std::endl;
                newPhasePresence = nPhaseOnly;

                //switch not depending on formulation
                //switch "Sw" to "xn20"
                priVars[switchIdx]
                        = volVars.moleFraction(nPhaseIdx, wCompIdx /*H2O*/);

                //switch all secondary components to mole fraction in non-wetting phase
                for (int compIdx=numMajorComponents; compIdx<numComponents; ++compIdx)
                    priVars[compIdx] = volVars.moleFraction(nPhaseIdx,compIdx);
            }
            //if saturation of non-wetting phase is smaller than 0 switch
            else if (volVars.saturation(nPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                //non-wetting phase has to disappear
                std::cout << "Non-wetting Phase disappears at vertex " << dofIdxGlobal
                            << ", coordinated: " << globalPos << ", Sn: "
                            << volVars.saturation(nPhaseIdx) << std::endl;
                newPhasePresence = wPhaseOnly;

                //switch "Sn" to "xwN2"
                priVars[switchIdx] = volVars.moleFraction(wPhaseIdx, nCompIdx /*N2*/);
            }
        }
        else if (phasePresence == nPhaseOnly)
        {
            Scalar xwmax = 1;
            Scalar sumxw = 0;
            //Calculate sum of mole fractions in the hypothetical wetting phase
            for (int compIdx = 0; compIdx < numComponents; compIdx++)
            {
                sumxw += volVars.moleFraction(wPhaseIdx, compIdx);
            }
            if (sumxw > xwmax)
                wouldSwitch = true;
            if (this->wasSwitched_[dofIdxGlobal])
                xwmax *=1.02;
            //wetting phase appears if sum is larger than one
            if (sumxw/*sum of mole fractions*/ > xwmax/*1*/)
            {
                std::cout << "Wetting Phase appears at vertex " << dofIdxGlobal
                        << ", coordinated: " << globalPos << ", sumxw: "
                        << sumxw << std::endl;
                newPhasePresence = bothPhases;

                //saturation of the wetting phase set to 0.0001 (if formulation TwoPFormulation::pnsw and vice versa)
                if (formulation == TwoPFormulation::pnsw)
                    priVars[switchIdx] = 0.0001;
                else if (formulation == TwoPFormulation::pwsn)
                    priVars[switchIdx] = 0.9999;

                //switch all secondary components back to wetting mole fraction
                for (int compIdx=numMajorComponents; compIdx<numComponents; ++compIdx)
                    priVars[compIdx] = volVars.moleFraction(wPhaseIdx,compIdx);
            }
        }
        else if (phasePresence == wPhaseOnly)
        {
            Scalar xnmax = 1;
            Scalar sumxn = 0;
            //Calculate sum of mole fractions in the hypothetical wetting phase
            for (int compIdx = 0; compIdx < numComponents; compIdx++)
            {
                sumxn += volVars.moleFraction(nPhaseIdx, compIdx);
            }
            if (sumxn > xnmax)
                wouldSwitch = true;
            if (this->wasSwitched_[dofIdxGlobal])
                xnmax *=1.02;
            //wetting phase appears if sum is larger than one
            if (sumxn > xnmax)
            {
                std::cout << "Non-wetting Phase appears at vertex " << dofIdxGlobal
                        << ", coordinated: " << globalPos << ", sumxn: "
                        << sumxn << std::endl;
                newPhasePresence = bothPhases;
                //saturation of the wetting phase set to 0.9999 (if formulation TwoPFormulation::pnsw and vice versa)
                if (formulation == TwoPFormulation::pnsw)
                    priVars[switchIdx] = 0.9999;
                else if (formulation == TwoPFormulation::pwsn)
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
