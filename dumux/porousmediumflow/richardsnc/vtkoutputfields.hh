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
 * \brief Adds vtk output fields specific to the Richards model.
 */
#ifndef DUMUX_RICHARDSNC_VTK_OUTPUT_FIELDS_HH
#define DUMUX_RICHARDSNC_VTK_OUTPUT_FIELDS_HH

#include <dumux/common/parameters.hh>

namespace Dumux {

/*!
 * \ingroup RichardsNCModel
 * \brief Adds vtk output fields specific to the Richards model.
 */
class RichardsNCVtkOutputFields
{
public:
    template <class VtkOutputModule>
    static void init(VtkOutputModule& vtk)
    {
        using VolumeVariables = typename VtkOutputModule::VolumeVariables;
        using FluidSystem = typename VolumeVariables::FluidSystem;

        vtk.addVolumeVariable([](const auto& v){ return v.saturation(VolumeVariables::liquidPhaseIdx); }, "S_w");
        vtk.addVolumeVariable([](const auto& v){ return v.saturation(VolumeVariables::gasPhaseIdx); }, "S_n");
        vtk.addVolumeVariable([](const auto& v){ return v.pressure(VolumeVariables::liquidPhaseIdx); }, "p_w");
        vtk.addVolumeVariable([](const auto& v){ return v.pressure(VolumeVariables::gasPhaseIdx); }, "p_n");
        vtk.addVolumeVariable([](const auto& v){ return v.capillaryPressure(); }, "pc");
        vtk.addVolumeVariable([](const auto& v){ return v.density(VolumeVariables::liquidPhaseIdx); }, "density");
        vtk.addVolumeVariable([](const auto& v){ return v.mobility(VolumeVariables::liquidPhaseIdx); }, "mobility");
        vtk.addVolumeVariable([](const auto& v){ return v.relativePermeability(VolumeVariables::liquidPhaseIdx); }, "kr");
        vtk.addVolumeVariable([](const auto& v){ return v.porosity(); }, "porosity");
        vtk.addVolumeVariable([](const auto& v){ return v.temperature(); }, "temperature");

        static const bool gravity = getParamFromGroup<bool>(vtk.paramGroup(), "Problem.EnableGravity");
        if(gravity)
            vtk.addVolumeVariable([](const auto& v){ return v.pressureHead(VolumeVariables::liquidPhaseIdx); }, "pressure head");
        vtk.addVolumeVariable([](const auto& v){ return v.waterContent(VolumeVariables::liquidPhaseIdx); },"water content");

        for (int compIdx = 0; compIdx < VolumeVariables::numComponents(); ++compIdx)
            vtk.addVolumeVariable([compIdx](const auto& v){ return v.moleFraction(VolumeVariables::liquidPhaseIdx, compIdx); },
                                  "x^" + FluidSystem::componentName(compIdx) + "_" + FluidSystem::phaseName(VolumeVariables::liquidPhaseIdx));

    }
};

} // end namespace Dumux

#endif
