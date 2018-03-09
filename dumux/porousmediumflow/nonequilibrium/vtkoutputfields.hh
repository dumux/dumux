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
 * \ingroup PorousmediumNonEquilibriumModel
 * \brief Adds vtk output fields specific to non-isothermal models
 */
#ifndef DUMUX_NONEQUILBRIUM_OUTPUT_FIELDS_HH
#define DUMUX_NONEQUILBRIUM_OUTPUT_FIELDS_HH

#include <dumux/common/properties.hh>

namespace Dumux
{

/*!
 * \ingroup PorousmediumNonEquilibriumModel
 * \brief Adds vtk output fields specific to non-isothermal models
 */
template<class EquilibriumVtkOutputFields, class FluidSystem, int numEnergyEqFluid, int numEnergyEqSolid>
class NonEquilibriumVtkOutputFields
{

public:
    template <class VtkOutputModule>
    static void init(VtkOutputModule& vtk)
    {
        EquilibriumVtkOutputFields::init(vtk);
         for (int i = 0; i < numEnergyEqFluid; ++i)
        vtk.addVolumeVariable( [i](const auto& v){ return v.temperature(i); }, "T_" + FluidSystem::phaseName(i) );
         for (int i = 0; i < numEnergyEqSolid; ++i)
        vtk.addVolumeVariable( [i](const auto& v){ return v.temperatureSolid(); }, "T_solid" );
        for (int i = 0; i < FluidSystem::numPhases; ++i){
        vtk.addVolumeVariable( [i](const auto& v){ return v.reynoldsNumber(i); }, "reynoldsNumber_" + FluidSystem::phaseName(i) );
        vtk.addVolumeVariable( [i](const auto& v){ return v.nusseltNumber(i); }, "nusseltNumber_" + FluidSystem::phaseName(i) );
        vtk.addVolumeVariable( [i](const auto& v){ return v.prandtlNumber(i); }, "prandtlNumber_" + FluidSystem::phaseName(i) );
        }
    }
};

} // end namespace Dumux

#endif
