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

#include <dumux/porousmediumflow/richards/vtkoutputfields.hh>
namespace Dumux
{

/*!
 * \ingroup RichardsNCModel
 * \brief Adds vtk output fields specific to the Richards model.
 */
template<class FluidSystem, class Indices>
class RichardsNCVtkOutputFields
{

public:
    template <class VtkOutputModule>
    static void init(VtkOutputModule& vtk)
    {
        vtk.addVolumeVariable([](const auto& v){ return v.saturation(Indices::wPhaseIdx); }, "Sw");
        vtk.addVolumeVariable([](const auto& v){ return v.saturation(Indices::nPhaseIdx); }, "Sn");
        vtk.addVolumeVariable([](const auto& v){ return v.pressure(Indices::wPhaseIdx); }, "pw");
        vtk.addVolumeVariable([](const auto& v){ return v.pressure(Indices::nPhaseIdx); }, "pn");
        vtk.addVolumeVariable([](const auto& v){ return v.capillaryPressure(); }, "pc");
        vtk.addVolumeVariable([](const auto& v){ return v.density(Indices::wPhaseIdx); }, "density");
        vtk.addVolumeVariable([](const auto& v){ return v.mobility(Indices::wPhaseIdx); }, "mobility");
        vtk.addVolumeVariable([](const auto& v){ return v.relativePermeability(Indices::wPhaseIdx); }, "kr");
        vtk.addVolumeVariable([](const auto& v){ return v.porosity(); }, "porosity");
        vtk.addVolumeVariable([](const auto& v){ return v.temperature(); }, "temperature");

        static const bool gravity = getParamFromGroup<bool>(vtk.paramGroup(), "Problem.EnableGravity");
        if(gravity)
            vtk.addVolumeVariable([](const auto& v){ return v.pressureHead(Indices::wPhaseIdx); }, "pressure head");
        vtk.addVolumeVariable([](const auto& v){ return v.waterContent(Indices::wPhaseIdx); },"water content");

        for (int k = 0; k < FluidSystem::numComponents; ++k)
            vtk.addVolumeVariable([k](const auto& v){ return v.moleFraction(Indices::wPhaseIdx, k); }, "x^" + FluidSystem::phaseName(Indices::wPhaseIdx) + "_" + FluidSystem::componentName(k));

    }
};

} // end namespace Dumux

#endif
