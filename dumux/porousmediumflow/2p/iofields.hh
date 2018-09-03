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
 * \brief Adds I/O fields specific to the two-phase model
 */
#ifndef DUMUX_TWOP_IO_FIELDS_HH
#define DUMUX_TWOP_IO_FIELDS_HH

#include <dune/common/deprecated.hh>

#include <dumux/io/name.hh>

#include "model.hh"

namespace Dumux {

/*!
 * \ingroup TwoPModel
 * \brief Adds I/O fields specific to the two-phase model
 */
template <TwoPFormulation priVarFormulation>
class TwoPIOFields
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        using VolumeVariables = typename VtkOutputModule::VolumeVariables;
        using FS = typename VolumeVariables::FluidSystem;

        out.addVolumeVariable([](const VolumeVariables& v){ return v.saturation(FS::phase0Idx); }, "S_"+FS::phaseName(FS::phase0Idx));
        out.addVolumeVariable([](const VolumeVariables& v){ return v.saturation(FS::phase1Idx); }, "S_"+FS::phaseName(FS::phase1Idx));
        out.addVolumeVariable([](const VolumeVariables& v){ return v.pressure(FS::phase0Idx); }, "p_"+FS::phaseName(FS::phase0Idx));
        out.addVolumeVariable([](const VolumeVariables& v){ return v.pressure(FS::phase1Idx); }, "p_"+FS::phaseName(FS::phase1Idx));
        out.addVolumeVariable([](const auto& v){ return v.capillaryPressure(); }, "pc");
        out.addVolumeVariable([](const auto& v){ return v.density(FS::phase0Idx); }, "rho_"+FS::phaseName(FS::phase0Idx));
        out.addVolumeVariable([](const auto& v){ return v.density(FS::phase1Idx); }, "rho_"+FS::phaseName(FS::phase1Idx));
        out.addVolumeVariable([](const auto& v){ return v.mobility(FS::phase0Idx); },"mob_"+ FS::phaseName(FS::phase0Idx));
        out.addVolumeVariable([](const auto& v){ return v.mobility(FS::phase1Idx); },"mob_"+ FS::phaseName(FS::phase1Idx));

        out.addVolumeVariable([](const auto& v){ return v.porosity(); }, "porosity");
    }

    template <class OutputModule>
    DUNE_DEPRECATED_MSG("use initOutputModule instead")
    static void init(OutputModule& out)
    {
        initOutputModule(out);
    }

    template <class FluidSystem = void, class SolidSystem = void>
    static std::string primaryVariableName(int pvIdx, int state = 0)
    {
        if (priVarFormulation == TwoPFormulation::p0s1)
            return pvIdx == 0 ? "p_w" : "S_n";
        else
            return pvIdx == 0 ? "p_n" : "S_w";
    }
};

} // end namespace Dumux

#endif
