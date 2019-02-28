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
 * \ingroup RichardsModel
 * \brief Adds I/O fields specific to the Richards model.
 */

#ifndef DUMUX_RICHARDS_IO_FIELDS_HH
#define DUMUX_RICHARDS_IO_FIELDS_HH

#include <dumux/common/parameters.hh>
#include <dumux/io/name.hh>

namespace Dumux {

/*!
 * \ingroup RichardsModel
 * \brief Adds I/O fields specific to the Richards model.
 */
template <bool enableWaterDiffusionInAir>
class RichardsIOFields
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        using VolumeVariables = typename OutputModule::VolumeVariables;
        using FS = typename VolumeVariables::FluidSystem;

        out.addVolumeVariable([](const auto& v){ return v.saturation(FS::liquidPhaseIdx); },
                              IOName::saturation<FS>(FS::liquidPhaseIdx));
        out.addVolumeVariable([](const auto& v){ return v.saturation(FS::gasPhaseIdx); },
                              IOName::saturation<FS>(FS::gasPhaseIdx));
        out.addVolumeVariable([](const auto& v){ return v.pressure(FS::liquidPhaseIdx); },
                              IOName::pressure<FS>(FS::liquidPhaseIdx));
        out.addVolumeVariable([](const auto& v){ return v.pressure(FS::gasPhaseIdx); },
                              IOName::pressure<FS>(FS::gasPhaseIdx));
        out.addVolumeVariable([](const auto& v){ return v.capillaryPressure(); },
                              IOName::capillaryPressure());
        out.addVolumeVariable([](const auto& v){ return v.density(FS::liquidPhaseIdx); },
                              IOName::density<FS>(FS::liquidPhaseIdx));
        out.addVolumeVariable([](const auto& v){ return v.mobility(FS::liquidPhaseIdx); },
                              IOName::mobility<FS>(FS::liquidPhaseIdx));
        out.addVolumeVariable([](const auto& v){ return v.relativePermeability(FS::liquidPhaseIdx); },
                              IOName::relativePermeability<FS>(FS::liquidPhaseIdx));
        out.addVolumeVariable([](const auto& v){ return v.porosity(); },
                              IOName::porosity());

        static const bool gravity = getParamFromGroup<bool>(out.paramGroup(), "Problem.EnableGravity");

        if(gravity)
            out.addVolumeVariable([](const auto& v){ return v.pressureHead(FS::liquidPhaseIdx); },
                                  IOName::pressureHead());
        if (enableWaterDiffusionInAir)
            out.addVolumeVariable([](const auto& v){ return v.moleFraction(FS::gasPhaseIdx, FS::liquidCompIdx); },
                                  IOName::moleFraction<FS>(FS::gasPhaseIdx, FS::liquidCompIdx));
        out.addVolumeVariable([](const auto& v){ return v.waterContent(FS::liquidPhaseIdx); },
                              IOName::waterContent());

        out.addVolumeVariable([](const auto& v){ return v.priVars().state(); },
                              IOName::phasePresence());
    }

    template<class ModelTraits, class FluidSystem, class SolidSystem = void>
    static std::string primaryVariableName(int pvIdx, int state)
    {
        using Indices = typename ModelTraits::Indices;

        if (state == Indices::gasPhaseOnly)
            return IOName::moleFraction<FluidSystem>(FluidSystem::gasPhaseIdx,
                                                     FluidSystem::liquidCompIdx);
        else
            return IOName::pressure<FluidSystem>(FluidSystem::liquidPhaseIdx);
    }
};

} // end namespace Dumux

#endif
