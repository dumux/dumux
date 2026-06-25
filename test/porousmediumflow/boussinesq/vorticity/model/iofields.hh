// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_BOUSSINESQ_VORTICITY_IO_FIELDS_HH
#define DUMUX_BOUSSINESQ_VORTICITY_IO_FIELDS_HH

#include <string>

namespace Dumux {

struct BoussinesqVorticityIOFields
{
    template<class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        using VV = typename OutputModule::VolumeVariables;
        using MT = typename VV::ModelTraits;

        constexpr int nPot  = MT::numPotentialEqs;
        constexpr int nComp = MT::numFluidComponents();

        if constexpr (nPot == 1)
            out.addVolumeVariable([](const auto& v){ return v.vectorPotential(0); }, "streamfunction");
        else
            for (int k = 0; k < nPot; ++k)
                out.addVolumeVariable(
                    [k](const auto& v){ return v.vectorPotential(k); },
                    "vectorPotential_" + std::to_string(k));

        if constexpr (nComp == 1)
            out.addVolumeVariable([](const auto& v){ return v.concentration(0); }, "concentration");
        else
            for (int i = 0; i < nComp; ++i)
                out.addVolumeVariable(
                    [i](const auto& v){ return v.concentration(i); },
                    "concentration_" + std::to_string(i));

        out.addVolumeVariable([](const auto& v){ return v.porosity(); }, "porosity");
    }

    template<class ModelTraits, class FS = void, class SolidSystem = void>
    static std::string primaryVariableName(int pvIdx, int /*state*/ = 0)
    {
        constexpr int nPot = ModelTraits::numPotentialEqs;
        if (pvIdx < nPot)
        {
            if (nPot == 1) return "streamfunction";
            return "vectorPotential_" + std::to_string(pvIdx);
        }
        const int i = pvIdx - nPot;
        if (ModelTraits::numFluidComponents() == 1) return "concentration";
        return "concentration_" + std::to_string(i);
    }
};

} // namespace Dumux

#endif
