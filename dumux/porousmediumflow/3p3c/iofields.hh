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
 * \brief Adds I/O fields specific to the three-phase three-component model
 */
#ifndef DUMUX_THREEPTHREEC_IO_FIELDS_HH
#define DUMUX_THREEPTHREEC_IO_FIELDS_HH

#include <dumux/io/name.hh>

namespace Dumux {

/*!
 * \ingroup ThreePThreeCModel
 * \brief Adds I/O fields specific to the three-phase three-component model
 */
template <class Indices>
class ThreePThreeCIOFields
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        using VolumeVariables = typename OutputModule::VolumeVariables;
        using FluidSystem = typename VolumeVariables::FluidSystem;

        // register standardized output fields
        for (int phaseIdx = 0; phaseIdx < VolumeVariables::numPhases(); ++phaseIdx)
        {
            out.addVolumeVariable( [phaseIdx](const auto& v){ return v.saturation(phaseIdx); },
                                   IOName::saturation<FluidSystem>(phaseIdx));
            out.addVolumeVariable( [phaseIdx](const auto& v){ return v.pressure(phaseIdx); },
                                   IOName::pressure<FluidSystem>(phaseIdx));
            out.addVolumeVariable( [phaseIdx](const auto& v){ return v.density(phaseIdx); },
                                   IOName::density<FluidSystem>(phaseIdx));

            for (int compIdx = 0; compIdx < VolumeVariables::numComponents(); ++compIdx)
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

    template <class OutputModule>
    DUNE_DEPRECATED_MSG("use initOutputModule instead")
    static void init(OutputModule& out)
    {
        initOutputModule(out);
    }

    template <class FluidSystem, class SolidSystem = void>
    static std::string primaryVariableName(int pvIdx, int state)
    {
        switch (state)
        {
            case Indices::threePhases:
            {
                const std::vector<std::string> s1 = {IOName::pressure<FluidSystem>(FluidSystem::gPhaseIdx),
                                                     IOName::saturation<FluidSystem>(FluidSystem::wPhaseIdx),
                                                     IOName::saturation<FluidSystem>(FluidSystem::nPhaseIdx)};
                return s1[pvIdx];
            }
            case Indices::wPhaseOnly:
            {
                const std::vector<std::string> s2 = {IOName::pressure<FluidSystem>(FluidSystem::gPhaseIdx),
                                                     IOName::moleFraction<FluidSystem>(FluidSystem::wPhaseIdx, FluidSystem::gCompIdx),
                                                     IOName::moleFraction<FluidSystem>(FluidSystem::wPhaseIdx, FluidSystem::nCompIdx)};
                return s2[pvIdx];
            }
            case Indices::gnPhaseOnly:
            {
                const std::vector<std::string> s3 = {IOName::pressure<FluidSystem>(FluidSystem::gPhaseIdx),
                                                     IOName::moleFraction<FluidSystem>(FluidSystem::gPhaseIdx, FluidSystem::wCompIdx),
                                                     IOName::saturation<FluidSystem>(FluidSystem::nPhaseIdx)};
                return s3[pvIdx];
            }
            case Indices::wnPhaseOnly:
            {
                const std::vector<std::string> s4 = {IOName::pressure<FluidSystem>(FluidSystem::gPhaseIdx),
                                                     IOName::moleFraction<FluidSystem>(FluidSystem::wPhaseIdx, FluidSystem::gCompIdx),
                                                     IOName::saturation<FluidSystem>(FluidSystem::nPhaseIdx)};
                return s4[pvIdx];
            }
            case Indices::gPhaseOnly:
            {
                const std::vector<std::string> s5 = {IOName::pressure<FluidSystem>(FluidSystem::gPhaseIdx),
                                                     IOName::moleFraction<FluidSystem>(FluidSystem::gPhaseIdx, FluidSystem::wCompIdx),
                                                     IOName::moleFraction<FluidSystem>(FluidSystem::gPhaseIdx, FluidSystem::nCompIdx)};
                return s5[pvIdx];
            }
            case Indices::wgPhaseOnly:
            {
                const std::vector<std::string> s6 = {IOName::pressure<FluidSystem>(FluidSystem::gPhaseIdx),
                                                     IOName::saturation<FluidSystem>(FluidSystem::wPhaseIdx),
                                                     IOName::moleFraction<FluidSystem>(FluidSystem::gPhaseIdx, FluidSystem::nCompIdx)};
                return s6[pvIdx];
            }
        }
    }
};

} // end namespace Dumux

#endif
