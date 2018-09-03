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
 * \brief Adds I/O fields specific to the OnePNC model
 */
#ifndef DUMUX_ONEPNC_IO_FIELDS_HH
#define DUMUX_ONEPNC_IO_FIELDS_HH

#include <string>

namespace Dumux {

/*!
 * \ingroup OnePNCModel
 * \brief Adds I/O fields specific to the OnePNC model
 */
template <bool useMoles>
class OnePNCIOFields
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        using VolumeVariables = typename OutputModule::VolumeVariables;
        using FluidSystem = typename VolumeVariables::FluidSystem;

        out.addVolumeVariable([](const auto& volVars){ return volVars.pressure(0); }, "p");
        out.addVolumeVariable([](const auto& volVars){ return volVars.density(0); }, "rho");
        out.addVolumeVariable([](const auto& volVars){ return volVars.viscosity(0); }, "mu");
        out.addVolumeVariable([](const auto& volVars){ return volVars.pressure(0) - 1e5; }, "delp");

        for (int i = 0; i < VolumeVariables::numComponents(); ++i)
           out.addVolumeVariable([i](const auto& volVars){ return volVars.moleFraction(0, i); },
                                     "x^" + std::string(FluidSystem::componentName(i)) + "_" + std::string(FluidSystem::phaseName(0)));

        for (int i = 0; i < VolumeVariables::numComponents(); ++i)
           out.addVolumeVariable([i](const auto& volVars){ return volVars.massFraction(0, i); },
                                     "X^" + std::string(FluidSystem::componentName(i))+ "_" + std::string(FluidSystem::phaseName(0)));
    }

    template <class OutputModule>
    DUNE_DEPRECATED_MSG("use initOutputModule instead")
    static void init(OutputModule& out)
    {
        initOutputModule(out);
    }

    template <class FluidSystem, class SolidSystem = void>
    static std::string primaryVariableName(int pvIdx, int state = 0)
    {
        const std::string xString = useMoles ? "x" : "X";
        if (pvIdx == 0)
            return "p";
        else
            return xString + "^" + FluidSystem::componentName(pvIdx)
                   + "_" + FluidSystem::phaseName(0);
    }
};

} // end namespace Dumux

#endif
