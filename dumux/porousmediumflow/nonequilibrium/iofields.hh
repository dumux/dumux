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
 * \ingroup NonEquilibriumModel
 * \brief Adds I/O fields specific to non-equilibrium models
 */

#ifndef DUMUX_NONEQUILBRIUM_OUTPUT_FIELDS_HH
#define DUMUX_NONEQUILBRIUM_OUTPUT_FIELDS_HH

#include <dumux/io/name.hh>

namespace Dumux {

template<class ModelTraits, class EquilibriumIOFields, bool enableThermalNonEquilibrium>
class NonEquilibriumIOFieldsImplementation;

template<class ModelTraits, class EquilibriumIOFields>
using NonEquilibriumIOFields =  NonEquilibriumIOFieldsImplementation<ModelTraits, EquilibriumIOFields, ModelTraits::enableThermalNonEquilibrium()>;
/*!
 * \ingroup NonEquilibriumModel
 * \brief Adds I/O fields specific to non-equilibrium models with chemical and thermal nonequilbirum or thermal non-equilibrium only
 */
template<class ModelTraits, class EquilibriumIOFields>
class NonEquilibriumIOFieldsImplementation<ModelTraits, EquilibriumIOFields, true>
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        using FluidSystem = typename OutputModule::VolumeVariables::FluidSystem;

        EquilibriumIOFields::initOutputModule(out);
        for (int i = 0; i < ModelTraits::numEnergyEqFluid(); ++i)
        {
            out.addVolumeVariable([i](const auto& v){ return v.temperatureFluid(i); },
            IOName::fluidTemperature<FluidSystem>(i));
        }

        out.addVolumeVariable([](const auto& v){ return v.temperatureSolid(); },
                              IOName::solidTemperature());

        for (int i = 0; i < ModelTraits::numFluidPhases(); ++i)
        {
            out.addVolumeVariable( [i](const auto& v){ return v.reynoldsNumber(i); }, "reynoldsNumber_" + FluidSystem::phaseName(i) );
            out.addVolumeVariable( [i](const auto& v){ return v.nusseltNumber(i); }, "nusseltNumber_" + FluidSystem::phaseName(i) );
            out.addVolumeVariable( [i](const auto& v){ return v.prandtlNumber(i); }, "prandtlNumber_" + FluidSystem::phaseName(i) );
        }
    }
};

template<class ModelTraits, class EquilibriumIOFields>
class NonEquilibriumIOFieldsImplementation<ModelTraits, EquilibriumIOFields, false>
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        using FluidSystem = typename OutputModule::VolumeVariables::FluidSystem;

        EquilibriumIOFields::initOutputModule(out);

        for (int i = 0; i < ModelTraits::numFluidPhases(); ++i)
        {
            out.addVolumeVariable( [i](const auto& v){ return v.reynoldsNumber(i); }, "reynoldsNumber_" + FluidSystem::phaseName(i) );
        }
    }
};

} // end namespace Dumux

#endif
