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
 * \brief Adds I/O fields specific to the three-phase three-component model.
 */

#ifndef DUMUX_THREEPTHREEC_IO_FIELDS_HH
#define DUMUX_THREEPTHREEC_IO_FIELDS_HH

#include <array>
#include <dumux/io/name.hh>

namespace Dumux {

/*!
 * \ingroup ThreePThreeCModel
 * \brief Adds I/O fields specific to the three-phase three-component model.
 */
class ThreePThreeCIOFields
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        using VolumeVariables = typename OutputModule::VolumeVariables;
        using FluidSystem = typename VolumeVariables::FluidSystem;

        // register standardized output fields
        for (int phaseIdx = 0; phaseIdx < VolumeVariables::numFluidPhases(); ++phaseIdx)
        {
            out.addVolumeVariable( [phaseIdx](const auto& v){ return v.saturation(phaseIdx); },
                                   IOName::saturation<FluidSystem>(phaseIdx));
            out.addVolumeVariable( [phaseIdx](const auto& v){ return v.pressure(phaseIdx); },
                                   IOName::pressure<FluidSystem>(phaseIdx));
            out.addVolumeVariable( [phaseIdx](const auto& v){ return v.density(phaseIdx); },
                                   IOName::density<FluidSystem>(phaseIdx));

            for (int compIdx = 0; compIdx < VolumeVariables::numFluidComponents(); ++compIdx)
                out.addVolumeVariable([phaseIdx, compIdx](const auto& v){ return v.moleFraction(phaseIdx, compIdx); },
                                      IOName::moleFraction<FluidSystem>(phaseIdx, compIdx));
        }

        out.addVolumeVariable( [](const auto& v){ return v.porosity(); },
                               IOName::porosity());
        out.addVolumeVariable( [](const auto& v){ return v.priVars().state(); },
                               IOName::phasePresence());
        out.addVolumeVariable( [](const auto& v){ return v.permeability(); },
                               IOName::permeability());
    }

    template <class ModelTraits, class FluidSystem, class SolidSystem = void>
    static std::string primaryVariableName(int pvIdx, int state)
    {
        using Indices = typename ModelTraits::Indices;
        static constexpr auto numEq = ModelTraits::numEq();
        using StringVec = std::array<std::string, numEq>;

        switch (state)
        {
            case Indices::threePhases:
            {
                static const StringVec s1 = {IOName::pressure<FluidSystem>(FluidSystem::gPhaseIdx),
                                             IOName::saturation<FluidSystem>(FluidSystem::wPhaseIdx),
                                             IOName::saturation<FluidSystem>(FluidSystem::nPhaseIdx)};
                return s1[pvIdx];
            }
            case Indices::wPhaseOnly:
            {
                static const StringVec s2 = {IOName::pressure<FluidSystem>(FluidSystem::gPhaseIdx),
                                             IOName::moleFraction<FluidSystem>(FluidSystem::wPhaseIdx, FluidSystem::gCompIdx),
                                             IOName::moleFraction<FluidSystem>(FluidSystem::wPhaseIdx, FluidSystem::nCompIdx)};
                return s2[pvIdx];
            }
            case Indices::gnPhaseOnly:
            {
                static const StringVec s3 = {IOName::pressure<FluidSystem>(FluidSystem::gPhaseIdx),
                                             IOName::moleFraction<FluidSystem>(FluidSystem::gPhaseIdx, FluidSystem::wCompIdx),
                                             IOName::saturation<FluidSystem>(FluidSystem::nPhaseIdx)};
                return s3[pvIdx];
            }
            case Indices::wnPhaseOnly:
            {
                static const StringVec s4 = {IOName::pressure<FluidSystem>(FluidSystem::gPhaseIdx),
                                             IOName::moleFraction<FluidSystem>(FluidSystem::wPhaseIdx, FluidSystem::gCompIdx),
                                             IOName::saturation<FluidSystem>(FluidSystem::nPhaseIdx)};
                return s4[pvIdx];
            }
            case Indices::gPhaseOnly:
            {
                static const StringVec s5 = {IOName::pressure<FluidSystem>(FluidSystem::gPhaseIdx),
                                             IOName::moleFraction<FluidSystem>(FluidSystem::gPhaseIdx, FluidSystem::wCompIdx),
                                             IOName::moleFraction<FluidSystem>(FluidSystem::gPhaseIdx, FluidSystem::nCompIdx)};
                return s5[pvIdx];
            }
            case Indices::wgPhaseOnly:
            {
                static const StringVec s6 = {IOName::pressure<FluidSystem>(FluidSystem::gPhaseIdx),
                                             IOName::saturation<FluidSystem>(FluidSystem::wPhaseIdx),
                                             IOName::moleFraction<FluidSystem>(FluidSystem::gPhaseIdx, FluidSystem::nCompIdx)};
                return s6[pvIdx];
            }
        }
    }
};

} // end namespace Dumux

#endif
