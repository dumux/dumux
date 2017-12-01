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
 * \brief Adds vtk output fields specific to the twop-nc-min model
 */
#ifndef DUMUX_TWOP_NC_MIN_VTK_OUTPUT_FIELDS_HH
#define DUMUX_TWOP_NC_MIN_VTK_OUTPUT_FIELDS_HH

#include <dumux/porousmediumflow/2pnc/implicit/vtkoutputfields.hh>

namespace Dumux
{

/*!
 * \ingroup TwoPNCMin, InputOutput
 * \brief Adds vtk output fields specific to the TwoPNCMin model
 */
template<class TypeTag>
class TwoPNCMinVtkOutputFields
{
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

    static constexpr int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);
    static constexpr int numSPhases = GET_PROP_VALUE(TypeTag, NumSPhases);
    static constexpr int numComponents = GET_PROP_VALUE(TypeTag, NumComponents);
    static constexpr int dim = GridView::dimension;

public:
    template <class VtkOutputModule>
    static void init(VtkOutputModule& vtk)
    {
        // use default fields from the 2pnc model
        TwoPNCVtkOutputFields<TypeTag>::init(vtk);

        //output additional to TwoPNC output:
        for (int i = 0; i < numSPhases; ++i)
        {
            vtk.addVolumeVariable([i](const VolumeVariables& v){ return v.precipitateVolumeFraction(numPhases + i); },"precipVolFrac_"+ FluidSystem::phaseName(numPhases + i));
        }

        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.permeability()[0][0]; },"Kxx");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.permeability()[1][1]; },"Kyy");
    }

};

} // end namespace Dumux

#endif
