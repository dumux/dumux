// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePModel
 * \brief Adds I/O fields specific to the one phase model.
 */

#ifndef DUMUX_ONEP_IO_FIELDS_HH
#define DUMUX_ONEP_IO_FIELDS_HH

#include <dumux/io/name.hh>

namespace Dumux {

/*!
 * \ingroup OnePModel
 * \brief Adds I/O fields specific to the one phase model.
 */
class OnePIOFields
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        out.addVolumeVariable([](const auto& volVars){ return volVars.pressure(); },
                              IOName::pressure());
    }

    template <class ModelTraits = void, class FluidSystem = void, class SolidSystem = void>
    static std::string primaryVariableName(int pvIdx = 0, int state = 0)
    {
        return IOName::pressure();
    }
};

} // end namespace Dumux

#endif
