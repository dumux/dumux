// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePNCModel
 * \brief Adds I/O fields specific to the OnePNC model.
 */

#ifndef DUMUX_ONEPNC_IO_FIELDS_HH
#define DUMUX_ONEPNC_IO_FIELDS_HH

#include <string>
#include <dumux/io/name.hh>

namespace Dumux {

/*!
 * \ingroup OnePNCModel
 * \brief Adds I/O fields specific to the OnePNC model.
 */
class OnePNCIOFields
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        using VolumeVariables = typename OutputModule::VolumeVariables;
        using FluidSystem = typename VolumeVariables::FluidSystem;

        out.addVolumeVariable([](const auto& volVars){ return volVars.pressure(0); },
                              IOName::pressure());
        out.addVolumeVariable([](const auto& volVars){ return volVars.density(0); },
                              IOName::density());
        out.addVolumeVariable([](const auto& volVars){ return volVars.viscosity(0); },
                              IOName::viscosity());

        for (int i = 0; i < VolumeVariables::numFluidComponents(); ++i)
           out.addVolumeVariable([i](const auto& volVars){ return volVars.moleFraction(0, i); },
                                     IOName::moleFraction<FluidSystem>(0, i));

        for (int i = 0; i < VolumeVariables::numFluidComponents(); ++i)
           out.addVolumeVariable([i](const auto& volVars){ return volVars.massFraction(0, i); },
                                     IOName::massFraction<FluidSystem>(0, i));
    }

    template <class ModelTraits, class FluidSystem, class SolidSystem = void>
    static std::string primaryVariableName(int pvIdx, int state = 0)
    {
        if (pvIdx == 0)
            return IOName::pressure();
        else if (ModelTraits::useMoles())
            return IOName::moleFraction<FluidSystem>(0, pvIdx);
        else
            return IOName::massFraction<FluidSystem>(0, pvIdx);
    }
};

} // end namespace Dumux

#endif
