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
 * \ingroup OnePNCModel
 * \brief Adds vtk output fields specific to the OnePNC model
 */
#ifndef DUMUX_ONEPNC_VTK_OUTPUT_FIELDS_HH
#define DUMUX_ONEPNC_VTK_OUTPUT_FIELDS_HH

namespace Dumux {

/*!
 * \ingroup OnePNCModel
 * \brief Adds vtk output fields specific to the OnePNC model
 */
class OnePNCVtkOutputFields
{
public:
    template <class VtkOutputModule>
    static void init(VtkOutputModule& vtk)
    {
        using VolumeVariables = typename VtkOutputModule::VolumeVariables;
        using FluidSystem = typename VolumeVariables::FluidSystem;
        using Indices = typename VolumeVariables::Indices;

        vtk.addVolumeVariable([](const auto& volVars){ return volVars.pressure(Indices::fluidSystemPhaseIdx); }, "pressure");
        vtk.addVolumeVariable([](const auto& volVars){ return volVars.density(Indices::fluidSystemPhaseIdx); }, "rho");
        vtk.addVolumeVariable([](const auto& volVars){ return volVars.viscosity(Indices::fluidSystemPhaseIdx); }, "mu");
        vtk.addVolumeVariable([](const auto& volVars){ return volVars.pressure(Indices::fluidSystemPhaseIdx) - 1e5; }, "delp");

        for (int i = 0; i < VolumeVariables::numComponents(); ++i)
           vtk.addVolumeVariable([i](const auto& volVars){ return volVars.moleFraction(Indices::fluidSystemPhaseIdx, i); },
                                     "x_" + std::string(FluidSystem::componentName(i)));

        for (int i = 0; i < VolumeVariables::numComponents(); ++i)
           vtk.addVolumeVariable([i](const auto& volVars){ return volVars.massFraction(Indices::fluidSystemPhaseIdx, i); },
                                     "X_" + std::string(FluidSystem::componentName(i)));
    }
};

} // end namespace Dumux

#endif
