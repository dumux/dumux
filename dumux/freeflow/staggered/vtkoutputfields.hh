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
#ifndef DUMUX_NAVIER_STOKES_VTK_OUTPUT_FIELDS_HH
#define DUMUX_NAVIER_STOKES_VTK_OUTPUT_FIELDS_HH

#include <dumux/common/basicproperties.hh>

namespace Dumux
{

/*!
 * \ingroup TwoP, InputOutput
 * \brief Adds vtk output fields specific to the twop model
 */
template<class TypeTag>
class NavierStokesVtkOutputFields
{
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
public:
    template <class VtkOutputModule>
    static void init(VtkOutputModule& vtk)
    {
        // vtk.addSecondaryVariable("Sw", [](const VolumeVariables& v){ return v.saturation(Indices::wPhaseIdx); });
        // vtk.addSecondaryVariable("Sn", [](const VolumeVariables& v){ return v.saturation(Indices::nPhaseIdx); });
        // vtk.addSecondaryVariable("pw", [](const VolumeVariables& v){ return v.pressure(Indices::wPhaseIdx); });
        // vtk.addSecondaryVariable("pn", [](const VolumeVariables& v){ return v.pressure(Indices::nPhaseIdx); });
        // vtk.addSecondaryVariable("p", [](const VolumeVariables& v){ return v.pressure(); });
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.pressure(); }, "p");
        // vtk.addSecondaryVariable("pc", [](const VolumeVariables& v){ return v.capillaryPressure(); });
        // vtk.addSecondaryVariable("rhoW", [](const VolumeVariables& v){ return v.density(Indices::wPhaseIdx); });
        // vtk.addSecondaryVariable("rhoN", [](const VolumeVariables& v){ return v.density(Indices::nPhaseIdx); });
        // vtk.addSecondaryVariable("mobW", [](const VolumeVariables& v){ return v.mobility(Indices::wPhaseIdx); });
        // vtk.addSecondaryVariable("mobN", [](const VolumeVariables& v){ return v.mobility(Indices::nPhaseIdx); });
        // vtk.addSecondaryVariable("temperature", [](const VolumeVariables& v){ return v.temperature(); });
        // vtk.addSecondaryVariable("porosity", [](const VolumeVariables& v){ return v.porosity(); });
    }
};

} // end namespace Dumux

#endif
