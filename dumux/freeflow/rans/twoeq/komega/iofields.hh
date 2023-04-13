// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup KOmegaModel
 * \copydoc Dumux::KOmegaIOFields
 */
#ifndef DUMUX_KOMEGA_IO_FIELDS_HH
#define DUMUX_KOMEGA_IO_FIELDS_HH

#include <dumux/freeflow/rans/iofields.hh>

namespace Dumux {

/*!
 * \ingroup KOmegaModel
 * \brief Adds I/O fields for the Reynolds-Averaged Navier-Stokes model
 */
struct KOmegaIOFields
{
    //! Initialize the KOmegaModel specific output fields.
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        RANSIOFields::initOutputModule(out);

        out.addVolumeVariable([](const auto& v){ return v.turbulentKineticEnergy(); }, "k");
        out.addVolumeVariable([](const auto& v){ return v.dissipation(); }, "omega");
    }

    //! return the names of the primary variables
    template <class ModelTraits, class FluidSystem>
    static std::string primaryVariableName(int pvIdx = 0, int state = 0)
    {
        if (pvIdx < ModelTraits::dim() + ModelTraits::numFluidComponents())
            return RANSIOFields::template primaryVariableName<ModelTraits, FluidSystem>(pvIdx, state);
        else if (pvIdx == ModelTraits::dim() + ModelTraits::numFluidComponents())
            return "k";
        else
            return "omega";
    }
};

} // end namespace Dumux

#endif
