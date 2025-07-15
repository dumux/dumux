// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesMassOnePIOFields
 */
#ifndef DUMUX_NAVIERSTOKES_MASS_1P_IO_FIELDS_HH
#define DUMUX_NAVIERSTOKES_MASS_1P_IO_FIELDS_HH

#include <dumux/common/parameters.hh>
#include <dumux/io/name.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Adds I/O fields for the Navier-Stokes model
 */
class NavierStokesMassOnePIOFields
{
public:
    //! Initialize the Navier-Stokes specific output fields.
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        out.addVolumeVariable([](const auto& v){ return v.pressure(); }, IOName::pressure());
        out.addVolumeVariable([](const auto& v){ return v.density(); }, IOName::density());
    }

    //! return the names of the primary variables
    template <class ModelTraits, class FluidSystem = void>
    static std::string primaryVariableName(int pvIdx = 0, int state = 0)
    {
        return IOName::pressure();
    }
};

} // end namespace Dumux

#endif
