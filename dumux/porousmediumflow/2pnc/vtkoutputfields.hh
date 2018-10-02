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

#include <dumux/porousmediumflow/2p/vtkoutputfields.hh>

namespace Dumux
{

/*!
 * \ingroup TwoPNCModel
 * \brief Adds vtk output fields specific to the TwoPNC model
 */
class TwoPNCVtkOutputFields
{
public:
   template <class VtkOutputModule>
    static void init(VtkOutputModule& vtk)
    {
        using VolumeVariables = typename VtkOutputModule::VolumeVariables;
        using FluidSystem = typename VolumeVariables::FluidSystem;

        // use default fields from the 2p model
        TwoPVtkOutputFields::init(vtk);

        //output additional to TwoP output:
        for (int i = 0; i < VolumeVariables::numPhases(); ++i)
            for (int j = 0; j < VolumeVariables::numComponents(); ++j)
                vtk.addVolumeVariable([i,j](const auto& v){ return v.moleFraction(i,j); },
                                    "x^"+ FluidSystem::componentName(j) + "_" + FluidSystem::phaseName(i));

        for (int i = 0; i < VolumeVariables::numPhases(); ++i)
            vtk.addVolumeVariable([i](const auto& v){ return v.molarDensity(i); },
                                    "rhoMolar_" + FluidSystem::phaseName(i));

        if (VolumeVariables::numComponents() < 3){
            for (int i = 0; i < VolumeVariables::numPhases(); ++i)
                for (int j = 0; j < VolumeVariables::numComponents(); ++j)
                    vtk.addVolumeVariable([i,j](const auto& v){ return v.massFraction(i,j); },
                                    "X^"+ FluidSystem::componentName(j) + "_" + FluidSystem::phaseName(i));
        }

        vtk.addVolumeVariable([](const auto& v){ return v.priVars().state(); }, "phase presence");
    }
};


} // end namespace Dumux

#endif
