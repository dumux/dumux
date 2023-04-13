// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PNMTwoPModel
 * \copydoc Dumux::PoreNetwork::TwoPIOFields
 */
#ifndef DUMUX_PNM_TWOP_IO_FIELDS_HH
#define DUMUX_PNM_TWOP_IO_FIELDS_HH

#include <dumux/porousmediumflow/2p/iofields.hh>
#include <dumux/porenetwork/common/iofields.hh>

namespace Dumux::PoreNetwork {

/*!
 * \ingroup PNMTwoPModel
 * \brief Adds output fields specific to the PNM 2p model
 */
class TwoPIOFields
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
         // use default fields from the 2p model
        Dumux::TwoPIOFields::initOutputModule(out);
        using VolumeVariables = typename OutputModule::VolumeVariables;
        using FS = typename VolumeVariables::FluidSystem;

        CommonIOFields::initOutputModule(out);

        out.addFluxVariable([](const auto& fluxVars, const auto& fluxVarsCache)
                              { return fluxVarsCache.pcEntry(); }, "pcEntry");

        out.addFluxVariable([](const auto& fluxVars, const auto& fluxVarsCache)
                              { return fluxVarsCache.transmissibility(FS::phase0Idx); }, "transmissibilityW");

        out.addFluxVariable([](const auto& fluxVars, const auto& fluxVarsCache)
                              { return fluxVarsCache.transmissibility(FS::phase1Idx); }, "transmissibilityN");

        auto volumeFluxW = [](const auto& fluxVars, const auto& fluxVarsCache)
        {
            auto upwindTerm = [](const auto& volVars) { return volVars.mobility(FS::phase0Idx); };
            using std::abs;
            return abs(fluxVars.advectiveFlux(FS::phase0Idx, upwindTerm));
        };
        out.addFluxVariable(volumeFluxW, "volumeFluxW");

        auto volumeFluxN = [](const auto& fluxVars, const auto& fluxVarsCache)
        {
            auto upwindTerm = [](const auto& volVars) { return volVars.mobility(FS::phase1Idx); };
            using std::abs;
            return abs(fluxVars.advectiveFlux(FS::phase1Idx, upwindTerm));
        };
        out.addFluxVariable(volumeFluxN, "volumeFluxN");
    }

    template <class ModelTraits, class FluidSystem, class SolidSystem = void>
    static std::string primaryVariableName(int pvIdx, int state = 0)
    {
        return TwoPIOFields::template primaryVariableName<ModelTraits, FluidSystem, SolidSystem>(pvIdx, state);
    }
};

} // end namespace Dumux::PoreNetwork

#endif
