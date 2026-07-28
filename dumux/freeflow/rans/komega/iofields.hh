// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \copydoc Dumux::KOmegaMassIOFields
 */
#ifndef DUMUX_RANS_KOMEGA_MASS_IO_FIELDS_HH
#define DUMUX_RANS_KOMEGA_MASS_IO_FIELDS_HH

#include <dumux/freeflow/navierstokes/mass/1p/iofields.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowModels
 * \brief Adds I/O fields for the two-equation k-omega turbulence closure fused onto the
 *        single-phase Navier-Stokes mass model.
 */
class KOmegaMassIOFields
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        NavierStokesMassOnePIOFields::initOutputModule(out);
        out.addVolumeVariable([](const auto& v){ return v.turbulentKineticEnergy(); }, "k");
        out.addVolumeVariable([](const auto& v){ return v.dissipation(); }, "omega");
        out.addVolumeVariable([](const auto& v){ return v.dynamicEddyViscosity(); }, "mu_t");
    }

    template <class ModelTraits, class FluidSystem = void>
    static std::string primaryVariableName(int pvIdx = 0, int state = 0)
    {
        if (pvIdx == ModelTraits::Indices::pressureIdx)
            return NavierStokesMassOnePIOFields::template primaryVariableName<ModelTraits, FluidSystem>(pvIdx, state);
        else if (pvIdx == ModelTraits::Indices::turbulentKineticEnergyIdx)
            return "k";
        else
            return "omega";
    }
};

} // end namespace Dumux

#endif
