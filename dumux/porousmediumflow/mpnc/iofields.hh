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
 * \ingroup MPNCModel
 * \brief Adds I/O fields specific to the twop model
 */
#ifndef DUMUX_MPNC_IO_FIELDS_HH
#define DUMUX_MPNC_IO_FIELDS_HH

namespace Dumux {

/*!
 * \ingroup MPNCModel
 * \brief Adds I/O fields specific to the twop model
 */
class MPNCIOFields
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        using VolumeVariables = typename OutputModule::VolumeVariables;
        using FluidSystem = typename VolumeVariables::FluidSystem;

        for (int i = 0; i < VolumeVariables::numPhases(); ++i)
        {
            out.addVolumeVariable([i](const auto& v){ return v.saturation(i); }, "S_"+ FluidSystem::phaseName(i));
            out.addVolumeVariable([i](const auto& v){ return v.pressure(i); }, "p_"+ FluidSystem::phaseName(i));
            out.addVolumeVariable([i](const auto& v){ return v.density(i); }, "rho_"+ FluidSystem::phaseName(i));
            out.addVolumeVariable([i](const auto& v){ return v.mobility(i); },"mob_"+ FluidSystem::phaseName(i));

            for (int j = 0; j < VolumeVariables::numComponents(); ++j)
                vtk.addVolumeVariable([i,j](const auto& v){ return v.moleFraction(i,j); },
                                      "x^"+ FluidSystem::componentName(j) + "_" + FluidSystem::phaseName(i));
        }
        for (int j = 0; j < VolumeVariables::numComponents(); ++j)
            vtk.addVolumeVariable([j](const auto& v){ return v.fugacity(j); },
                                  "fugacity^"+ FluidSystem::componentName(j));
    }
};

} // end namespace Dumux

#endif
