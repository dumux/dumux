// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoroElastic
 * \brief Adds I/O fields specific to the poro-elastic model
 */
#ifndef DUMUX_POROELASTIC_IO_FIELDS_HH
#define DUMUX_POROELASTIC_IO_FIELDS_HH

#include <dumux/io/name.hh>

namespace Dumux {

/*!
 * \ingroup PoroElastic
 * \brief Adds I/O fields specific to the poro-elastic model
 */
class PoroElasticIOFields
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        out.addVolumeVariable([](const auto& volVars){ return volVars.displacement(); },
                              IOName::displacement());
        out.addVolumeVariable([](const auto& volVars){ return volVars.porosity(); },
                              IOName::porosity());
    }

};

} // end namespace Dumux

#endif
