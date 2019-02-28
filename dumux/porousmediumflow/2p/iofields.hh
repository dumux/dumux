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
 * \ingroup TwoPModel
 * \brief Adds I/O fields specific to the two-phase model.
 */

#ifndef DUMUX_TWOP_IO_FIELDS_HH
#define DUMUX_TWOP_IO_FIELDS_HH

#include <dumux/io/name.hh>

#include "model.hh"

namespace Dumux {

/*!
 * \ingroup TwoPModel
 * \brief Adds I/O fields specific to the two-phase model.
 */
class TwoPIOFields
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        using VolumeVariables = typename OutputModule::VolumeVariables;
        using FS = typename VolumeVariables::FluidSystem;

        for (int phaseIdx = 0; phaseIdx < FS::numPhases; ++phaseIdx)
        {
            out.addVolumeVariable([phaseIdx](const VolumeVariables& v){ return v.saturation(phaseIdx); },
                                  IOName::saturation<FS>(phaseIdx));
            out.addVolumeVariable([phaseIdx](const VolumeVariables& v){ return v.pressure(phaseIdx); },
                                  IOName::pressure<FS>(phaseIdx));
            out.addVolumeVariable([phaseIdx](const auto& v){ return v.density(phaseIdx); },
                                  IOName::density<FS>(phaseIdx));
            out.addVolumeVariable([phaseIdx](const auto& v){ return v.mobility(phaseIdx); },
                                  IOName::mobility<FS>(phaseIdx));
        }

        out.addVolumeVariable([](const auto& v){ return v.capillaryPressure(); },
                              IOName::capillaryPressure());
        out.addVolumeVariable([](const auto& v){ return v.porosity(); },
                              IOName::porosity());
    }

    template <class ModelTraits, class FluidSystem, class SolidSystem = void>
    static std::string primaryVariableName(int pvIdx, int state = 0)
    {
        if (ModelTraits::priVarFormulation() == TwoPFormulation::p0s1)
            return pvIdx == 0 ? IOName::pressure<FluidSystem>(FluidSystem::phase0Idx)
                              : IOName::saturation<FluidSystem>(FluidSystem::phase1Idx);
        else
            return pvIdx == 0 ? IOName::pressure<FluidSystem>(FluidSystem::phase1Idx)
                              : IOName::saturation<FluidSystem>(FluidSystem::phase0Idx);
    }
};

} // end namespace Dumux

#endif
