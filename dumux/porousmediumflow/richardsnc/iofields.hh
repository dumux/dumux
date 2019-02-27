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
 * \ingroup RichardsNCModel
 * \brief Adds I/O fields specific to the Richards model.
 */

#ifndef DUMUX_RICHARDSNC_IO_FIELDS_HH
#define DUMUX_RICHARDSNC_IO_FIELDS_HH

#include <dumux/common/parameters.hh>
#include <dumux/io/name.hh>

namespace Dumux {

/*!
 * \ingroup RichardsNCModel
 * \brief Adds I/O fields specific to the Richards model.
 */
class RichardsNCIOFields
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        using VolumeVariables = typename OutputModule::VolumeVariables;
        using FS = typename VolumeVariables::FluidSystem;

        out.addVolumeVariable([](const auto& v){ return v.saturation(VolumeVariables::liquidPhaseIdx); },
                              IOName::saturation() + "_" + IOName::liquidPhase());
        out.addVolumeVariable([](const auto& v){ return v.saturation(VolumeVariables::gasPhaseIdx); },
                              IOName::saturation() + "_" + IOName::gaseousPhase());
        out.addVolumeVariable([](const auto& v){ return v.pressure(VolumeVariables::liquidPhaseIdx); },
                              IOName::pressure() + "_" + IOName::liquidPhase());
        out.addVolumeVariable([](const auto& v){ return v.pressure(VolumeVariables::gasPhaseIdx); },
                              IOName::pressure() + "_" + IOName::gaseousPhase());
        out.addVolumeVariable([](const auto& v){ return v.capillaryPressure(); },
                              IOName::capillaryPressure());
        out.addVolumeVariable([](const auto& v){ return v.density(FS::liquidPhaseIdx); },
                              IOName::density() + "_" + IOName::liquidPhase());
        out.addVolumeVariable([](const auto& v){ return v.mobility(FS::liquidPhaseIdx); },
                              IOName::mobility() + "_" + IOName::liquidPhase());
        out.addVolumeVariable([](const auto& v){ return v.relativePermeability(VolumeVariables::liquidPhaseIdx); },
                              IOName::relativePermeability() + "_" + IOName::liquidPhase());
        out.addVolumeVariable([](const auto& v){ return v.porosity(); },
                              IOName::porosity());
        out.addVolumeVariable([](const auto& v){ return v.temperature(); },
                              IOName::temperature());

        static const bool gravity = getParamFromGroup<bool>(out.paramGroup(), "Problem.EnableGravity");
        if (gravity)
            out.addVolumeVariable([](const auto& v){ return v.pressureHead(VolumeVariables::liquidPhaseIdx); },
                                  IOName::pressureHead());
        out.addVolumeVariable([](const auto& v){ return v.waterContent(VolumeVariables::liquidPhaseIdx); },
                              IOName::waterContent());

        for (int compIdx = 0; compIdx < VolumeVariables::numFluidComponents(); ++compIdx)
            out.addVolumeVariable([compIdx](const auto& v){ return v.moleFraction(VolumeVariables::liquidPhaseIdx, compIdx); },
                                  IOName::moleFraction<FS>(VolumeVariables::liquidPhaseIdx, compIdx));

    }

    template <class ModelTraits, class FluidSystem, class SolidSystem = void>
    static std::string primaryVariableName(int pvIdx, int state = 0)
    {
        if (pvIdx == 0)
            return IOName::pressure<FluidSystem>(0);
        else
            return ModelTraits::useMoles() ? IOName::moleFraction<FluidSystem>(0, pvIdx)
                                           : IOName::massFraction<FluidSystem>(0, pvIdx);
    }
};

} // end namespace Dumux

#endif
