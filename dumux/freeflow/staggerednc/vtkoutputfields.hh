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
 * \brief Adds vtk output fields specific to the NavierStokesNC model
 */
#ifndef DUMUX_NAVIER_STOKES_NC_VTK_OUTPUT_FIELDS_HH
#define DUMUX_NAVIER_STOKES_NC_VTK_OUTPUT_FIELDS_HH

#include <dumux/common/basicproperties.hh>

namespace Dumux
{

/*!
 * \ingroup TwoP, InputOutput
 * \brief Adds vtk output fields specific to the NavierStokesNC model
 */
template<class TypeTag>
class NavierStokesNCVtkOutputFields
{
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

    static constexpr int numComponents = GET_PROP_VALUE(TypeTag, NumComponents);
    static constexpr int phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);

public:
    template <class VtkOutputModule>
    static void init(VtkOutputModule& vtk)
    {
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.pressure(); }, "p");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.molarDensity(); }, "rhoMolar");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.density(); }, "rho");

        for (int j = 0; j < numComponents; ++j)
        {
            vtk.addVolumeVariable([j](const VolumeVariables& v){ return v.massFraction(phaseIdx,j); }, "X^" + FluidSystem::componentName(j) + "_" + FluidSystem::phaseName(phaseIdx));
            vtk.addVolumeVariable([j](const VolumeVariables& v){ return v.moleFraction(phaseIdx,j); }, "x^" + FluidSystem::componentName(j) + "_" + FluidSystem::phaseName(phaseIdx));
        }

        if(GET_PROP_VALUE(TypeTag, EnableEnergyBalance))
            vtk.addVolumeVariable( [](const VolumeVariables& v){ return v.temperature(); },"temperature");
    }
};

} // end namespace Dumux

#endif
