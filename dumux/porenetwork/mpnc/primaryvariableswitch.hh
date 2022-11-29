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
 * \ingroup PNMMPNCModel
 * \brief The primary variable switch for the PNM MPNC model.
 */

#ifndef DUMUX_MPNC_PRIMARY_VARIABLE_SWITCH_HH
#define DUMUX_MPNC_PRIMARY_VARIABLE_SWITCH_HH

#include <iostream>
#include <dumux/porousmediumflow/compositional/primaryvariableswitch.hh>

namespace Dumux::PoreNetwork {

/*!
 * \ingroup PNMMPNCModel
 * \brief The primary variable switch controlling the phase presence state variable.
 */
class MPNCPrimaryVariableSwitch
: public Dumux::PrimaryVariableSwitch<MPNCPrimaryVariableSwitch>
{
    using ParentType = Dumux::PrimaryVariableSwitch<MPNCPrimaryVariableSwitch>;
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
        using Indices = typename VolumeVariables::Indices;

        static constexpr auto numPhases = VolumeVariables::numFluidPhases();
        static_assert(numPhases == 2, "PV switch only implemented for two phases");

        static const Scalar thresholdSw = getParam<Scalar>("MPNC.ThresholdSw", 0.2);

        const auto wPhaseIdx = volVars.fluidState().wettingPhase();
        const auto nPhaseIdx = 1 - wPhaseIdx;

        // evaluate primary variable switch
        bool wouldSwitch = false;
        const int phase = priVars.state();
        int newPhase = phase;

        //check if a primary variable switch is necessary
        if (newPhase == VolumeVariables::pwIsPrimaryVariableIndex)
        {
            if (volVars.saturation(wPhaseIdx) <= thresholdSw)
            {
                wouldSwitch = true;
                if (this->verbosity() > 1)
                    std::cout << "Wetting phase (" << FluidSystem::phaseName(wPhaseIdx) << ")"
                              << " is below threshhold at dof " << dofIdxGlobal
                              << ", coordinates: " << globalPos
                              << ", S_" << FluidSystem::phaseName(wPhaseIdx) << ": " << volVars.saturation(wPhaseIdx)
                              << std::endl;

                newPhase = VolumeVariables::pnIsPrimaryVariableIndex;
                priVars[Indices::p0Idx] = volVars.pressure(nPhaseIdx);
            }
        }
        else
        {
            if (volVars.saturation(wPhaseIdx) > thresholdSw)
            {
                wouldSwitch = true;
                if (this->verbosity() > 1)
                    std::cout << "Wetting phase (" << FluidSystem::phaseName(wPhaseIdx) << ")"
                              << " is above threshhold at dof " << dofIdxGlobal
                              << ", coordinates: " << globalPos
                              << ", S_" << FluidSystem::phaseName(wPhaseIdx) << ": " << volVars.saturation(wPhaseIdx)
                              << std::endl;

                newPhase = VolumeVariables::pwIsPrimaryVariableIndex;
                priVars[Indices::p0Idx] = volVars.pressure(wPhaseIdx);
            }
        }

        priVars.setState(newPhase);
        this->wasSwitched_[dofIdxGlobal] = wouldSwitch;
        return phase != newPhase;
    }
};

} // end namespace Dumux::PoreNetwork

#endif
