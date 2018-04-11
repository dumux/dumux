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
 * \ingroup TwoPModel
 * \brief Adds vtk output fields specific to the two-phase model
 */
#ifndef DUMUX_TWOP_VTK_OUTPUT_FIELDS_HH
#define DUMUX_TWOP_VTK_OUTPUT_FIELDS_HH

namespace Dumux {

/*!
 * \ingroup TwoPModel
 * \brief Adds vtk output fields specific to the two-phase model
 */
class TwoPVtkOutputFields
{
public:
    template <class VtkOutputModule>
    static void init(VtkOutputModule& vtk)
    {
        using VolumeVariables = typename VtkOutputModule::VolumeVariables;
        using FS = typename VolumeVariables::FluidSystem;

        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.porosity(); }, "porosity");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.capillaryPressure(); }, "pc");

        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.saturation(FS::phase0Idx); }, "Sw");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.pressure(FS::phase0Idx); }, "pw");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.density(FS::phase0Idx); }, "rhoW");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.mobility(FS::phase0Idx); }, "mobW");

        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.saturation(FS::phase1Idx); }, "Sn");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.pressure(FS::phase1Idx); }, "pn");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.density(FS::phase1Idx); }, "rhoN");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.mobility(FS::phase1Idx); }, "mobN");
    }
};

} // end namespace Dumux

#endif
