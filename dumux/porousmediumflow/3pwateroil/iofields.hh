// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ThreePWaterOilModel
 * \brief Adds I/O fields specific to the three-phase three-component model.
 */

#ifndef DUMUX_3P2CNI_IO_FIELDS_HH
#define DUMUX_3P2CNI_IO_FIELDS_HH

#include <dumux/io/name.hh>
#include <dumux/porousmediumflow/3p3c/iofields.hh>

namespace Dumux {

/*!
 * \ingroup ThreePWaterOilModel
 * \brief Adds I/O fields specific to the three-phase three-component model.
 */
class ThreePWaterOilIOFields
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
            out.addVolumeVariable([phaseIdx](const auto& v){ return v.mobility(phaseIdx); },
                                  IOName::mobility<FluidSystem>(phaseIdx));
            out.addVolumeVariable([phaseIdx](const auto& v){ return v.viscosity(phaseIdx); },
                                  IOName::viscosity<FluidSystem>(phaseIdx));

            for (int compIdx = 0; compIdx < VolumeVariables::numFluidComponents(); ++compIdx)
                out.addVolumeVariable([phaseIdx, compIdx](const auto& v){ return v.moleFraction(phaseIdx, compIdx); },
                                      IOName::moleFraction<FluidSystem>(phaseIdx, compIdx));
        }

        out.addVolumeVariable([](const auto& v){ return v.porosity(); },
                              IOName::porosity());
        out.addVolumeVariable([](const auto& v){ return v.priVars().state(); },
                              IOName::phasePresence());
        out.addVolumeVariable([](const auto& v){ return v.permeability(); },
                              IOName::permeability());
    }

    template <class ModelTraits, class FluidSystem, class SolidSystem = void>
    static std::string primaryVariableName(int pvIdx, int state)
    {
        using Indices = typename ModelTraits::Indices;
        static constexpr auto numEq = ModelTraits::numEq();
        using StringVec = std::array<std::string, numEq>;

        switch (state)
        {
            case Indices::threePhases:
            {
                static const StringVec s1 = {IOName::pressure<FluidSystem>(FluidSystem::gPhaseIdx),
                                             IOName::saturation<FluidSystem>(FluidSystem::wPhaseIdx),
                                             IOName::saturation<FluidSystem>(FluidSystem::nPhaseIdx)};
                return s1[pvIdx];
            }
            case Indices::wPhaseOnly:
            {
                static const StringVec s2 = {IOName::pressure<FluidSystem>(FluidSystem::wPhaseIdx),
                                             IOName::temperature(),
                                             IOName::moleFraction<FluidSystem>(FluidSystem::wPhaseIdx, FluidSystem::nCompIdx)};
                return s2[pvIdx];
            }
            case Indices::gnPhaseOnly:
            {
                static const StringVec s3 = {IOName::pressure<FluidSystem>(FluidSystem::gPhaseIdx),
                                             IOName::saturation<FluidSystem>(FluidSystem::nPhaseIdx),
                                             IOName::moleFraction<FluidSystem>(FluidSystem::nPhaseIdx, FluidSystem::wCompIdx)};
                return s3[pvIdx];
            }
            case Indices::wnPhaseOnly:
            {
                static const StringVec s4 = {IOName::pressure<FluidSystem>(FluidSystem::wPhaseIdx),
                                             IOName::temperature(),
                                             IOName::saturation<FluidSystem>(FluidSystem::nPhaseIdx)};
                return s4[pvIdx];
            }
            case Indices::gPhaseOnly:
            {
                static const StringVec s5 = {IOName::pressure<FluidSystem>(FluidSystem::gPhaseIdx),
                                             IOName::temperature(),
                                             IOName::moleFraction<FluidSystem>(FluidSystem::gPhaseIdx, FluidSystem::nCompIdx)};
                return s5[pvIdx];
            }
            case Indices::wgPhaseOnly:
            {
                static const StringVec s6 = {IOName::pressure<FluidSystem>(FluidSystem::gPhaseIdx),
                                             IOName::saturation<FluidSystem>(FluidSystem::wPhaseIdx),
                                             IOName::moleFraction<FluidSystem>(FluidSystem::gPhaseIdx, FluidSystem::nCompIdx)};
                return s6[pvIdx];
            }
        }
    }
};

} // end namespace Dumux

#endif
