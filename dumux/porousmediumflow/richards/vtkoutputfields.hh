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
 * \brief Adds vtk output fields specific to the Richards model.
 */
#ifndef DUMUX_RICHARDS_VTK_OUTPUT_FIELDS_HH
#define DUMUX_RICHARDS_VTK_OUTPUT_FIELDS_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

namespace Dumux {

/*!
 * \ingroup RichardsModel
 * \brief Adds vtk output fields specific to the Richards model.
 */
template<class TypeTag>
class RichardsVtkOutputFields
{
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);

    enum {
        nPhaseIdx = Indices::nPhaseIdx,
        wPhaseIdx = Indices::wPhaseIdx
    };

    static constexpr bool enableWaterDiffusionInAir = GET_PROP_VALUE(TypeTag, EnableWaterDiffusionInAir);

public:
    template <class VtkOutputModule>
    static void init(VtkOutputModule& vtk)
    {
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.saturation(Indices::wPhaseIdx); }, "Sw");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.saturation(Indices::nPhaseIdx); }, "Sn");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.pressure(Indices::wPhaseIdx); }, "pw");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.pressure(Indices::nPhaseIdx); }, "pn");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.capillaryPressure(); }, "pc");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.density(Indices::wPhaseIdx); }, "density");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.mobility(Indices::wPhaseIdx); }, "mobility");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.relativePermeability(wPhaseIdx); }, "kr");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.porosity(); }, "porosity");

        static const bool gravity = getParamFromGroup<bool>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "Problem.EnableGravity");
        if(gravity)
            vtk.addVolumeVariable([](const VolumeVariables& v){ return v.pressureHead(wPhaseIdx); }, "pressure head");
        if (enableWaterDiffusionInAir)
            vtk.addVolumeVariable([](const VolumeVariables& v){ return v.moleFraction(1, 0); }, "x^w_air");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.waterContent(wPhaseIdx); },"water content");

         vtk.addVolumeVariable([](const VolumeVariables& v){ return v.priVars().state(); }, "phasePresence");

    }
};

} // end namespace Dumux

#endif
