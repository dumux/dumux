// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterModels
 * \brief Add I/O fields specific to shallow water
 */
#ifndef DUMUX_FREEFLOW_SHALLOW_WATER_IO_FIELDS_HH
#define DUMUX_FREEFLOW_SHALLOW_WATER_IO_FIELDS_HH

#include <dumux/io/name.hh>
#include <string>

namespace Dumux {

/*!
 * \ingroup ShallowWaterModels
 * \brief Adds vtk output fields for the shallow water model
 */
class ShallowWaterIOFields
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        using VolumeVariables = typename OutputModule::VolumeVariables;

        out.addVolumeVariable([](const VolumeVariables& v){ return v.waterDepth(); }, "waterDepth");
        out.addVolumeVariable([](const VolumeVariables& v){ return v.velocity(0); }, "velocityX");
        out.addVolumeVariable([](const VolumeVariables& v){ return v.velocity(1); }, "velocityY");
        out.addVolumeVariable([](const VolumeVariables& v){ return v.bedSurface(); }, "bedSurface");
        out.addVolumeVariable([](const VolumeVariables& v){ return v.bedSurface() + v.waterDepth(); }, "freeSurface");
    }

    /* Restart is limited for shallow water models since only primary
    * variables are regarded. Important parameters like the bedSurface,
    * fricition,.. are missing so far.
    */
    template <class ModelTraits>
    static std::string primaryVariableName(int pvIdx)
    {
        std::string name;

        switch(pvIdx){
            case 0 : name = "waterDepth";
            case 1 : name = "velocityX";
            case 2 : name = "velocityY";
        }

        return name;
    }

};

} // end namespace Dumux

#endif
