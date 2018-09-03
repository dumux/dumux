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
 * \ingroup RichardsModel
 * \brief Adds I/O fields specific to the Richards model.
 */
#ifndef DUMUX_RICHARDS_IO_FIELDS_HH
#define DUMUX_RICHARDS_IO_FIELDS_HH

#include <dumux/common/parameters.hh>

namespace Dumux {

/*!
 * \ingroup RichardsModel
 * \brief Adds I/O fields specific to the Richards model.
 */
template<bool enableWaterDiffusionInAir, class Indices>
class RichardsIOFields
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        using VolumeVariables = typename VtkOutputModule::VolumeVariables;
        using FS = typename VolumeVariables::FluidSystem;

        out.addVolumeVariable([](const auto& v){ return v.saturation(FS::liquidPhaseIdx); }, "S_"+FS::phaseName(FS::liquidPhaseIdx));
        out.addVolumeVariable([](const auto& v){ return v.saturation(FS::gasPhaseIdx); }, "S_"+FS::phaseName(FS::gasPhaseIdx));
        out.addVolumeVariable([](const auto& v){ return v.pressure(FS::liquidPhaseIdx); }, "p_"+FS::phaseName(FS::liquidPhaseIdx));
        out.addVolumeVariable([](const auto& v){ return v.pressure(FS::gasPhaseIdx); }, "p_"+FS::phaseName(FS::gasPhaseIdx));
        out.addVolumeVariable([](const auto& v){ return v.capillaryPressure(); }, "pc");
        out.addVolumeVariable([](const auto& v){ return v.density(FS::liquidPhaseIdx); }, "rho_"+FS::phaseName(FS::liquidPhaseIdx));
        out.addVolumeVariable([](const auto& v){ return v.mobility(FS::liquidPhaseIdx); }, "mob_"+FS::phaseName(FS::liquidPhaseIdx));
        out.addVolumeVariable([](const auto& v){ return v.relativePermeability(FS::liquidPhaseIdx); }, "kr");
        out.addVolumeVariable([](const auto& v){ return v.porosity(); }, "porosity");

        static const bool gravity = getParamFromGroup<bool>(out.paramGroup(), "Problem.EnableGravity");

        if(gravity)
            out.addVolumeVariable([](const auto& v){ return v.pressureHead(FS::liquidPhaseIdx); }, "pressure head");
        if (enableWaterDiffusionInAir)
            out.addVolumeVariable([](const auto& v){ return v.moleFraction(1, 0); }, "x^"+FS::componentName(FS::gasCompIdx)+"_"+FS::phaseName(FS::liquidPhaseIdx));
        out.addVolumeVariable([](const auto& v){ return v.waterContent(FS::liquidPhaseIdx); },"water content");

        out.addVolumeVariable([](const auto& v){ return v.priVars().state(); }, "phase presence");
    }

    template <class OutputModule>
    DUNE_DEPRECATED_MSG("use initOutputModule instead")
    static void init(OutputModule& out)
    {
        initOutputModule(out);
    }

    template<class FluidSystem = void, class SolidSystem = void>
    static std::string primaryVariableName(int pvIdx, int state)
    {
        if (state == Indices::gasPhaseOnly)
            return "x^w_n";
        else
            return "p_w";
    }
};

} // end namespace Dumux

#endif
