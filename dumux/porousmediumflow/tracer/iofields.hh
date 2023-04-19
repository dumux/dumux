// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TracerModel
 * \brief Adds I/O fields specific to the tracer model.
 */

#ifndef DUMUX_TRACER_IO_FIELDS_HH
#define DUMUX_TRACER_IO_FIELDS_HH

#include <string>

#include <dumux/io/name.hh>

namespace Dumux {

/*!
 * \ingroup TracerModel
 * \brief Adds I/O fields specific to the tracer model.
 */
class TracerIOFields
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        using VolumeVariables = typename OutputModule::VolumeVariables;
        using FluidSystem = typename VolumeVariables::FluidSystem;

        // register standardized out output fields
        for (int compIdx = 0; compIdx < VolumeVariables::numFluidComponents(); ++compIdx)
        {
            out.addVolumeVariable([compIdx](const auto& v){ return v.moleFraction(0, compIdx); },
                                  "x^" + FluidSystem::componentName(compIdx));
            out.addVolumeVariable([compIdx](const auto& v){ return v.massFraction(0, compIdx); },
                                  "X^" + FluidSystem::componentName(compIdx));
        }
        out.addVolumeVariable( [](const auto& v){ return v.density(); }, IOName::density());
    }

    template <class ModelTraits, class FluidSystem, class SolidSystem = void>
    static std::string primaryVariableName(int pvIdx, int state = 0)
    {
        const std::string xString = ModelTraits::useMoles() ? "x" : "X";
        return xString + "^" + FluidSystem::componentName(pvIdx);
    }
};

} // end namespace Dumux

#endif
