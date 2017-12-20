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
 * \brief Adds vtk output fields specific to non-isothermal models
 */
#ifndef DUMUX_NONEQUILBRIUM_OUTPUT_FIELDS_HH
#define DUMUX_NONEQUILBRIUM_OUTPUT_FIELDS_HH

#include <dumux/common/properties.hh>

namespace Dumux
{

/*!
 * \ingroup NonEquilibrium, InputOutput
 * \brief Adds vtk output fields specific to non-isothermal models
 */
template<class TypeTag>
class NonEquilibriumVtkOutputFields
{
    using EquilibriumVtkOutputFields = typename GET_PROP_TYPE(TypeTag, EquilibriumVtkOutputFields);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

    static constexpr int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);
    static constexpr int numEnergyEqFluid = GET_PROP_VALUE(TypeTag, NumEnergyEqFluid);
    static constexpr int numEnergyEqSolid = GET_PROP_VALUE(TypeTag, NumEnergyEqSolid);
public:
    template <class VtkOutputModule>
    static void init(VtkOutputModule& vtk)
    {
        EquilibriumVtkOutputFields::init(vtk);
         for (int i = 0; i < numEnergyEqFluid; ++i)
        vtk.addVolumeVariable( [i](const VolumeVariables& v){ return v.temperature(i); }, "T_" + FluidSystem::phaseName(i) );
         for (int i = 0; i < numEnergyEqSolid; ++i)
        vtk.addVolumeVariable( [i](const VolumeVariables& v){ return v.temperatureSolid(); }, "T_solid" );
        for (int i = 0; i < numPhases; ++i){
        vtk.addVolumeVariable( [i](const VolumeVariables& v){ return v.reynoldsNumber(i); }, "reynoldsNumber_" + FluidSystem::phaseName(i) );
        vtk.addVolumeVariable( [i](const VolumeVariables& v){ return v.nusseltNumber(i); }, "nusseltNumber_" + FluidSystem::phaseName(i) );
        vtk.addVolumeVariable( [i](const VolumeVariables& v){ return v.prandtlNumber(i); }, "prandtlNumber_" + FluidSystem::phaseName(i) );
        }
    }
};

} // end namespace Dumux

#endif
