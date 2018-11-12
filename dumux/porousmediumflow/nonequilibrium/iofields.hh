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
 * \brief Adds I/O fields specific to non-isothermal models
 */
#ifndef DUMUX_NONEQUILBRIUM_OUTPUT_FIELDS_HH
#define DUMUX_NONEQUILBRIUM_OUTPUT_FIELDS_HH

#include <dumux/io/name.hh>

namespace Dumux {

/*!
 * \ingroup PorousmediumNonEquilibriumModel
 * \brief Adds I/O fields specific to non-isothermal models
 */
template<class ModelTraits, class EquilibriumIOFields>
class NonEquilibriumIOFields
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        using FluidSystem = typename OutputModule::VolumeVariables::FluidSystem;

        EquilibriumIOFields::initOutputModule(out);
        for (int i = 0; i < ModelTraits::numEnergyEqFluid(); ++i)
            out.addVolumeVariable([i](const auto& v){ return v.temperatureFluid(i); },
                                  IOName::fluidTemperature<FluidSystem>(i));
        out.addVolumeVariable([](const auto& v){ return v.temperatureSolid(); },
                              IOName::solidTemperature());
        for (int i = 0; i < ModelTraits::numPhases(); ++i){
            out.addVolumeVariable( [i](const auto& v){ return v.reynoldsNumber(i); }, "reynoldsNumber_" + FluidSystem::phaseName(i) );
            out.addVolumeVariable( [i](const auto& v){ return v.nusseltNumber(i); }, "nusseltNumber_" + FluidSystem::phaseName(i) );
            out.addVolumeVariable( [i](const auto& v){ return v.prandtlNumber(i); }, "prandtlNumber_" + FluidSystem::phaseName(i) );
        }
    }

    template <class OutputModule>
    DUNE_DEPRECATED_MSG("use initOutputModule instead")
    static void init(OutputModule& out)
    {
        initOutputModule(out);
    }
};

} // end namespace Dumux

#endif
