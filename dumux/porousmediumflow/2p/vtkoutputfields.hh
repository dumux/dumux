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

#ifdef DUMUXFENICS // gas phase = wetting phase
        const int wPhaseIdx = FS::phase1Idx; // gas (set in spatial params)
        const int nPhaseIdx = FS::phase0Idx; // liquid
#else
        const int wPhaseIdx = FS::phase0Idx;
        const int nPhaseIdx = FS::phase1Idx;
#endif

        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.porosity(); }, "porosity");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.capillaryPressure(); }, "pc");

        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.saturation(wPhaseIdx); }, "S_w");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.pressure(wPhaseIdx); }, "p_w");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.density(wPhaseIdx); }, "rho_w");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.mobility(wPhaseIdx); }, "mob_w");

        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.saturation(nPhaseIdx); }, "S_n");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.pressure(nPhaseIdx); }, "p_n");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.density(nPhaseIdx); }, "rho_n");
        vtk.addVolumeVariable([](const VolumeVariables& v){ return v.mobility(nPhaseIdx); }, "mob_n");
    }
};

} // end namespace Dumux

#endif
