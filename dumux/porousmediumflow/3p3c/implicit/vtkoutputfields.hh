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
 * \brief Adds vtk output fields specific to the twop model
 */
#ifndef DUMUX_THREEPTHREEC_VTK_OUTPUT_FIELDS_HH
#define DUMUX_THREEPTHREEC_VTK_OUTPUT_FIELDS_HH

#include <dumux/common/properties.hh>

namespace Dumux
{

/*!
 * \ingroup ThreePThreeC, InputOutput
 * \brief Adds vtk output fields specific to the three-phase three-component model
 */
template<class TypeTag>
class ThreePThreeCVtkOutputFields
{
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
public:
    template <class VtkOutputModule>
    static void init(VtkOutputModule& vtk)
    {
        // register standardized vtk output fields
        vtk.addVolumeVariable( [](const VolumeVariables& v){ return v.saturation(Indices::wPhaseIdx); }, "Sw");
        vtk.addVolumeVariable( [](const VolumeVariables& v){ return v.saturation(Indices::nPhaseIdx); },"Sn");
        vtk.addVolumeVariable( [](const VolumeVariables& v){ return v.saturation(Indices::gPhaseIdx); },"Sg");
        vtk.addVolumeVariable( [](const VolumeVariables& v){ return v.pressure(Indices::wPhaseIdx); },"pw");
        vtk.addVolumeVariable( [](const VolumeVariables& v){ return v.pressure(Indices::nPhaseIdx); },"pn");
        vtk.addVolumeVariable( [](const VolumeVariables& v){ return v.pressure(Indices::gPhaseIdx); },"pg");
        vtk.addVolumeVariable( [](const VolumeVariables& v){ return v.density(Indices::wPhaseIdx); },"rhow");
        vtk.addVolumeVariable( [](const VolumeVariables& v){ return v.density(Indices::nPhaseIdx); },"rhon");
        vtk.addVolumeVariable( [](const VolumeVariables& v){ return v.density(Indices::gPhaseIdx); },"rhog");
        vtk.addVolumeVariable( [](const VolumeVariables& v){ return v.moleFraction(Indices::wPhaseIdx, Indices::wCompIdx); },"x^H2O_w");
        vtk.addVolumeVariable( [](const VolumeVariables& v){ return v.moleFraction(Indices::wPhaseIdx, Indices::nCompIdx); },"x^mesitylene_w");
        vtk.addVolumeVariable( [](const VolumeVariables& v){ return v.moleFraction(Indices::wPhaseIdx, Indices::gCompIdx); },"x^Air_w");
        vtk.addVolumeVariable( [](const VolumeVariables& v){ return v.moleFraction(Indices::nPhaseIdx, Indices::wCompIdx); },"x^H2O_n");
        vtk.addVolumeVariable( [](const VolumeVariables& v){ return v.moleFraction(Indices::nPhaseIdx, Indices::nCompIdx); },"x^mesitylene_n");
        vtk.addVolumeVariable( [](const VolumeVariables& v){ return v.moleFraction(Indices::nPhaseIdx, Indices::gCompIdx); },"x^Air_n");
        vtk.addVolumeVariable( [](const VolumeVariables& v){ return v.moleFraction(Indices::gPhaseIdx, Indices::wCompIdx); },"x^H2O_g");
        vtk.addVolumeVariable( [](const VolumeVariables& v){ return v.moleFraction(Indices::gPhaseIdx, Indices::nCompIdx); },"x^mesitylene_g");
        vtk.addVolumeVariable( [](const VolumeVariables& v){ return v.moleFraction(Indices::gPhaseIdx, Indices::gCompIdx); },"x^Air_g");
        vtk.addVolumeVariable( [](const VolumeVariables& v){ return v.porosity(); },"porosity");
        vtk.addVolumeVariable( [](const VolumeVariables& v){ return v.temperature(); },"temperature");
        vtk.addVolumeVariable( [](const VolumeVariables& v){ return v.priVars().state(); },"phase presence");
        vtk.addVolumeVariable( [](const VolumeVariables& v){ return v.permeability(); },"permeability");
    }
};

} // end namespace Dumux

#endif
