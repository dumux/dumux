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
#ifndef DUMUX_THREEP_VTK_OUTPUT_FIELDS_HH
#define DUMUX_THREEP_VTK_OUTPUT_FIELDS_HH

#include <dumux/implicit/properties.hh>

namespace Dumux
{

/*!
 * \ingroup ThreeP, InputOutput
 * \brief Adds vtk output fields specific to the twop model
 */
template<class TypeTag>
class ThreePVtkOutputFields
{
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
public:
    template <class VtkOutputModule>
    static void init(VtkOutputModule& vtk)
    {
        // register standardized vtk output fields
        vtk.addSecondaryVariable("sw", [](const VolumeVariables& v){ return v.saturation(Indices::wPhaseIdx); });
        vtk.addSecondaryVariable("sn", [](const VolumeVariables& v){ return v.saturation(Indices::nPhaseIdx); });
        vtk.addSecondaryVariable("sg", [](const VolumeVariables& v){ return v.saturation(Indices::gPhaseIdx); });
        vtk.addSecondaryVariable("pw", [](const VolumeVariables& v){ return v.pressure(Indices::wPhaseIdx); });
        vtk.addSecondaryVariable("pn", [](const VolumeVariables& v){ return v.pressure(Indices::nPhaseIdx); });
        vtk.addSecondaryVariable("pg", [](const VolumeVariables& v){ return v.pressure(Indices::gPhaseIdx); });
        vtk.addSecondaryVariable("rhow", [](const VolumeVariables& v){ return v.density(Indices::wPhaseIdx); });
        vtk.addSecondaryVariable("rhon", [](const VolumeVariables& v){ return v.density(Indices::nPhaseIdx); });
        vtk.addSecondaryVariable("rhog", [](const VolumeVariables& v){ return v.density(Indices::gPhaseIdx); });
        vtk.addSecondaryVariable("porosity", [](const VolumeVariables& v){ return v.porosity(); });
        vtk.addSecondaryVariable("permeability", [](const VolumeVariables& v){ return v.permeability(); });
        vtk.addSecondaryVariable("temperature", [](const VolumeVariables& v){ return v.temperature(); });
    }
};

} // end namespace Dumux

#endif
