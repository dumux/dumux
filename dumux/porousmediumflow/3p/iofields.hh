// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup ThreePModel
 * \brief Adds I/O fields specific to the three-phase model.
 */

#ifndef DUMUX_THREEP_IO_FIELDS_HH
#define DUMUX_THREEP_IO_FIELDS_HH

#include <dumux/io/name.hh>

namespace Dumux {

/*!
 * \ingroup ThreePModel
 * \brief Adds I/O fields specific to the three-phase model.
 */
class ThreePIOFields
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        using VolumeVariables = typename OutputModule::VolumeVariables;
        using FluidSystem = typename VolumeVariables::FluidSystem;

        // register standardized output fields
        for (int phaseIdx = 0; phaseIdx < VolumeVariables::numFluidPhases(); ++phaseIdx)
        {
            out.addVolumeVariable([phaseIdx](const auto& v){ return v.saturation(phaseIdx); },
                                  IOName::saturation<FluidSystem>(phaseIdx));
            out.addVolumeVariable([phaseIdx](const auto& v){ return v.pressure(phaseIdx); },
                                  IOName::pressure<FluidSystem>(phaseIdx));
            out.addVolumeVariable([phaseIdx](const auto& v){ return v.density(phaseIdx); },
                                  IOName::density<FluidSystem>(phaseIdx));
        }

        out.addVolumeVariable( [](const auto& v){ return v.porosity(); },
                               IOName::porosity());
        out.addVolumeVariable( [](const auto& v){ return v.permeability(); },
                               IOName::permeability());
    }

    template <class ModelTraits, class FluidSystem, class SolidSystem = void>
    static std::string primaryVariableName(int pvIdx, int state = 0)
    {
        switch (pvIdx)
        {
            case 0: return IOName::pressure<FluidSystem>(FluidSystem::gPhaseIdx);
            case 1: return IOName::saturation<FluidSystem>(FluidSystem::wPhaseIdx);
            default: return IOName::saturation<FluidSystem>(FluidSystem::nPhaseIdx);
        }
    }
};

} // end namespace Dumux

#endif
