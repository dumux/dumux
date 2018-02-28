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
 * \ingroup TwoPNCModel
 * \brief Adds vtk output fields specific to the twop-nc model
 */
#ifndef DUMUX_TWOP_NC_VTK_OUTPUT_FIELDS_HH
#define DUMUX_TWOP_NC_VTK_OUTPUT_FIELDS_HH

#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/2p/vtkoutputfields.hh>

namespace Dumux
{

/*!
 * \ingroup TwoPNCModel
 * \brief Adds vtk output fields specific to the TwoPNC model
 */
template<class TypeTag>
class TwoPNCVtkOutputFields
{
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

    static constexpr int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);
    static constexpr int numComponents = GET_PROP_VALUE(TypeTag, NumComponents);

public:
    template <class VtkOutputModule>
    static void init(VtkOutputModule& vtk)
    {
        // use default fields from the 2p model
        TwoPVtkOutputFields<Indices>::init(vtk);

        //output additional to TwoP output:
        for (int i = 0; i < numPhases; ++i)
            for (int j = 0; j < numComponents; ++j)
                vtk.addVolumeVariable([i,j](const VolumeVariables& v){ return v.moleFraction(i,j); },"x_"+ FluidSystem::phaseName(i) + "^" + FluidSystem::componentName(j));

        for (int j = 0; j < numComponents; ++j)
            vtk.addVolumeVariable([j](const VolumeVariables& v){ return v.molarity(Indices::wPhaseIdx,j); },"m_"+ FluidSystem::phaseName(Indices::wPhaseIdx) + "^" + FluidSystem::componentName(j));

        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.priVars().state(); }, "phasePresence");
    }
};

} // end namespace Dumux

#endif
