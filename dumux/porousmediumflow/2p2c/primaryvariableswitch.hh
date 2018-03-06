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
 * \ingroup TwoPTwoCModel
 * \brief The primary variable switch for the 2p2c model
 */
#ifndef DUMUX_2P2C_PRIMARY_VARIABLE_SWITCH_HH
#define DUMUX_2P2C_PRIMARY_VARIABLE_SWITCH_HH

#include <dumux/porousmediumflow/compositional/primaryvariableswitch.hh>
#include "indices.hh" // for TwoPTwoCFormulation

namespace Dumux
{
/*!
 * \ingroup TwoPTwoCModel
 * \brief The primary variable switch controlling the phase presence state variable
 */
template<class TypeTag>
class TwoPTwoCPrimaryVariableSwitch
: public PrimaryVariableSwitch<typename GET_PROP_TYPE(TypeTag, FVGridGeometry), TwoPTwoCPrimaryVariableSwitch<TypeTag>>
{
    using ParentType = PrimaryVariableSwitch<typename GET_PROP_TYPE(TypeTag, FVGridGeometry), TwoPTwoCPrimaryVariableSwitch<TypeTag>>;
    friend ParentType;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using GlobalPosition = Dune::FieldVector<Scalar, GridView::dimensionworld>;

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    enum {
        switchIdx = Indices::switchIdx,

        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        wCompIdx = Indices::wCompIdx,
        nCompIdx = Indices::nCompIdx,

        wPhaseOnly = Indices::wPhaseOnly,
        nPhaseOnly = Indices::nPhaseOnly,
        bothPhases = Indices::bothPhases,

        pwsn = TwoPTwoCFormulation::pwsn,
        pnsw = TwoPTwoCFormulation::pnsw,
        formulation = GET_PROP_VALUE(TypeTag, Formulation)
    };

    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

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

        // check if a primary var switch is necessary
        if (phasePresence == nPhaseOnly)
        {
            // calculate mole fraction in the hypothetic wetting phase
            Scalar xww = volVars.moleFraction(wPhaseIdx, wCompIdx);
            Scalar xwn = volVars.moleFraction(wPhaseIdx, nCompIdx);

            Scalar xwMax = 1.0;
            if (xww + xwn > xwMax)
                wouldSwitch = true;
            if (this->wasSwitched_[dofIdxGlobal])
                xwMax *= 1.02;

            // if the sum of the mole fractions is larger than
            // 100%, wetting phase appears
            if (xww + xwn > xwMax)
            {
                // wetting phase appears
                std::cout << "wetting phase appears at vertex " << dofIdxGlobal
                          << ", coordinates: " << globalPos << ", xww + xwn: "
                          << xww + xwn << std::endl;
                newPhasePresence = bothPhases;
                if (formulation == pnsw)
                    priVars[switchIdx] = 0.0;
                else if (formulation == pwsn)
                    priVars[switchIdx] = 1.0;
            }
        }
        else if (phasePresence == wPhaseOnly)
        {
            // calculate fractions of the partial pressures in the
            // hypothetic nonwetting phase
            Scalar xnw = volVars.moleFraction(nPhaseIdx, wCompIdx);
            Scalar xnn = volVars.moleFraction(nPhaseIdx, nCompIdx);

            Scalar xgMax = 1.0;
            if (xnw + xnn > xgMax)
                wouldSwitch = true;
            if (this->wasSwitched_[dofIdxGlobal])
                xgMax *= 1.02;

            // if the sum of the mole fractions is larger than
            // 100%, nonwetting phase appears
            if (xnw + xnn > xgMax)
            {
                // nonwetting phase appears
                std::cout << "nonwetting phase appears at vertex " << dofIdxGlobal
                          << ", coordinates: " << globalPos << ", xnw + xnn: "
                          << xnw + xnn << std::endl;
                newPhasePresence = bothPhases;
                if (formulation == pnsw)
                    priVars[switchIdx] = 0.999;
                else if (formulation == pwsn)
                    priVars[switchIdx] = 0.001;
            }
        }
        else if (phasePresence == bothPhases)
        {
            Scalar Smin = 0.0;
            if (this->wasSwitched_[dofIdxGlobal])
                Smin = -0.01;

            if (volVars.saturation(nPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // nonwetting phase disappears
                std::cout << "Nonwetting phase disappears at vertex " << dofIdxGlobal
                          << ", coordinates: " << globalPos << ", sn: "
                          << volVars.saturation(nPhaseIdx) << std::endl;
                newPhasePresence = wPhaseOnly;

                if(useMoles) // mole-fraction formulation
                {
                    priVars[switchIdx]
                        = volVars.moleFraction(wPhaseIdx, nCompIdx);
                }
                else // mass-fraction formulation
                {
                    priVars[switchIdx]
                        = volVars.massFraction(wPhaseIdx, nCompIdx);
                }
            }
            else if (volVars.saturation(wPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // wetting phase disappears
                std::cout << "Wetting phase disappears at vertex " << dofIdxGlobal
                          << ", coordinates: " << globalPos << ", sw: "
                          << volVars.saturation(wPhaseIdx) << std::endl;
                newPhasePresence = nPhaseOnly;

                if(useMoles) // mole-fraction formulation
                {
                    priVars[switchIdx]
                        = volVars.moleFraction(nPhaseIdx, wCompIdx);
                }
                else // mass-fraction formulation
                {
                    priVars[switchIdx]
                        = volVars.massFraction(nPhaseIdx, wCompIdx);
                }
            }
        }

        priVars.setState(newPhasePresence);
        this->wasSwitched_[dofIdxGlobal] = wouldSwitch;
        return phasePresence != newPhasePresence;
    }
};

} // end namespace Dumux

#endif
