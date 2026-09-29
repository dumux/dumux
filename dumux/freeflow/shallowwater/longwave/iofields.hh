// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup LongWaveModel
 * \copydoc Dumux::LongWaveIOFields
 */
#ifndef DUMUX_FREEFLOW_SHALLOWWATER_LONGWAVE_IOFIELDS_HH
#define DUMUX_FREEFLOW_SHALLOWWATER_LONGWAVE_IOFIELDS_HH

#include <string>

namespace Dumux {

/*!
 * \ingroup LongWaveModel
 * \brief Adds output fields for the long-wave approximations of the shallow water equations
 */
class LongWaveIOFields
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        using VolumeVariables = typename OutputModule::VolumeVariables;

        out.addVolumeVariable([](const VolumeVariables& v){ return v.waterDepth(); }, "waterDepth");
        out.addVolumeVariable([](const VolumeVariables& v){ return v.bedSurface(); }, "bedSurface");
        out.addVolumeVariable([](const VolumeVariables& v){ return v.freeSurface(); }, "freeSurface");
    }

    template <class ModelTraits>
    static std::string primaryVariableName(int pvIdx = 0)
    { return "waterDepth"; }
};

} // end namespace Dumux

#endif
