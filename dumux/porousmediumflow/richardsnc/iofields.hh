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
 * \ingroup RichardsNCModel
 * \brief Adds I/O fields specific to the Richards model.
 */
#ifndef DUMUX_RICHARDSNC_IO_FIELDS_HH
#define DUMUX_RICHARDSNC_IO_FIELDS_HH

#include <dumux/common/parameters.hh>

namespace Dumux {

/*!
 * \ingroup RichardsNCModel
 * \brief Adds I/O fields specific to the Richards model.
 */
template <bool useMoles>
class RichardsNCIOFields
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        using VolumeVariables = typename OutputModule::VolumeVariables;
        using FS = typename VolumeVariables::FluidSystem;

        out.addVolumeVariable([](const auto& v){ return v.saturation(VolumeVariables::liquidPhaseIdx); }, "S_"+FS::phaseName(VolumeVariables::liquidPhaseIdx));
        out.addVolumeVariable([](const auto& v){ return v.saturation(VolumeVariables::gasPhaseIdx); }, "S_gas");
        out.addVolumeVariable([](const auto& v){ return v.pressure(VolumeVariables::liquidPhaseIdx); }, "p_"+FS::phaseName(VolumeVariables::liquidPhaseIdx));
        out.addVolumeVariable([](const auto& v){ return v.pressure(VolumeVariables::gasPhaseIdx); }, "p_gas");
        out.addVolumeVariable([](const auto& v){ return v.capillaryPressure(); }, "pc");
        out.addVolumeVariable([](const auto& v){ return v.density(FS::liquidPhaseIdx); }, "rho_"+FS::phaseName(FS::liquidPhaseIdx));
        out.addVolumeVariable([](const auto& v){ return v.mobility(FS::liquidPhaseIdx); }, "mob_"+FS::phaseName(FS::liquidPhaseIdx));
        out.addVolumeVariable([](const auto& v){ return v.relativePermeability(VolumeVariables::liquidPhaseIdx); }, "kr");
        out.addVolumeVariable([](const auto& v){ return v.porosity(); }, "porosity");
        out.addVolumeVariable([](const auto& v){ return v.temperature(); }, "T");

        static const bool gravity = getParamFromGroup<bool>(out.paramGroup(), "Problem.EnableGravity");
        if(gravity)
            out.addVolumeVariable([](const auto& v){ return v.pressureHead(VolumeVariables::liquidPhaseIdx); }, "pressure head");
        out.addVolumeVariable([](const auto& v){ return v.waterContent(VolumeVariables::liquidPhaseIdx); },"water content");

        for (int compIdx = 0; compIdx < VolumeVariables::numComponents(); ++compIdx)
            out.addVolumeVariable([compIdx](const auto& v){ return v.moleFraction(VolumeVariables::liquidPhaseIdx, compIdx); },
                                  "x^" + FS::componentName(compIdx) + "_" + FS::phaseName(VolumeVariables::liquidPhaseIdx));

    }

    template <class OutputModule>
    DUNE_DEPRECATED_MSG("use initOutputModule instead")
    static void init(OutputModule& out)
    {
        initOutputModule(out);
    }

    template <class FluidSystem, class SolidSystem = void>
    static std::string primaryVariableName(int pvIdx, int state = 0)
    {
        const std::string xString = useMoles ? "x" : "X";
        if (pvIdx == 0)
            return "p_" + FluidSystem::phaseName(0);
        else
            return xString + "^" + FluidSystem::componentName(pvIdx)
                   + "_" + FluidSystem::phaseName(0);
    }
};

} // end namespace Dumux

#endif
