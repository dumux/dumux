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
#include <dumux/porousmediumflow/2p/formulation.hh>

namespace Dumux
{
/*!
 * \ingroup TwoPTwoCModel
 * \brief The primary variable switch controlling the phase presence state variable
 */
class TwoPTwoCPrimaryVariableSwitch : public PrimaryVariableSwitch<TwoPTwoCPrimaryVariableSwitch>
{
    using ParentType = PrimaryVariableSwitch<TwoPTwoCPrimaryVariableSwitch>;
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
        static constexpr int wPhaseIdx = FluidSystem::wPhaseIdx;
        static constexpr int nPhaseIdx = FluidSystem::nPhaseIdx;
        static constexpr int wCompIdx = FluidSystem::wCompIdx;
        static constexpr int nCompIdx = FluidSystem::nCompIdx;

        static constexpr bool useMoles = VolumeVariables::useMoles();
        static constexpr auto formulation = VolumeVariables::priVarFormulation();

        using Indices = typename VolumeVariables::Indices;
        static constexpr int switchIdx = Indices::switchIdx;

        // evaluate primary variable switch
        bool wouldSwitch = false;
        int phasePresence = priVars.state();
        int newPhasePresence = phasePresence;

        // check if a primary var switch is necessary
        if (phasePresence == Indices::nPhaseOnly)
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
                newPhasePresence = Indices::bothPhases;
                if (formulation == TwoPFormulation::pnsw)
                    priVars[switchIdx] = 0.0001;
                else if (formulation == TwoPFormulation::pwsn)
                    priVars[switchIdx] = 0.9999;
            }
        }
        else if (phasePresence == Indices::wPhaseOnly)
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
                newPhasePresence = Indices::bothPhases;
                if (formulation == TwoPFormulation::pnsw)
                    priVars[switchIdx] = 0.9999;
                else if (formulation == TwoPFormulation::pwsn)
                    priVars[switchIdx] = 0.0001;
            }
        }
        else if (phasePresence == Indices::bothPhases)
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
                newPhasePresence = Indices::wPhaseOnly;

                if(useMoles) // mole-fraction formulation
                    priVars[switchIdx] = volVars.moleFraction(wPhaseIdx, nCompIdx);
                else // mass-fraction formulation
                    priVars[switchIdx] = volVars.massFraction(wPhaseIdx, nCompIdx);
            }
            else if (volVars.saturation(wPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // wetting phase disappears
                std::cout << "Wetting phase disappears at vertex " << dofIdxGlobal
                          << ", coordinates: " << globalPos << ", sw: "
                          << volVars.saturation(wPhaseIdx) << std::endl;
                newPhasePresence = Indices::nPhaseOnly;

                if(useMoles) // mole-fraction formulation
                    priVars[switchIdx] = volVars.moleFraction(nPhaseIdx, wCompIdx);
                else // mass-fraction formulation
                    priVars[switchIdx] = volVars.massFraction(nPhaseIdx, wCompIdx);
            }
        }

        priVars.setState(newPhasePresence);
        this->wasSwitched_[dofIdxGlobal] = wouldSwitch;
        return phasePresence != newPhasePresence;
    }
};

} // end namespace Dumux

#endif
