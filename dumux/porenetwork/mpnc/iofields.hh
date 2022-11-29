// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup PoreNetworkMNCPModel
 * \copydoc Dumux::PNMMPNCIOFields
 */
#ifndef DUMUX_PNM_MP_NC_IO_FIELDS_HH
#define DUMUX_PNM_MP_NC_IO_FIELDS_HH

#include <dumux/porousmediumflow/mpnc/iofields.hh>
#include <dumux/porenetwork/common/iofields.hh>

namespace Dumux::PoreNetwork {

/*!
 * \ingroup PoreNetworkOnePNCModel
 * \brief Adds output fields specific to the PNM MPNC model
 */
class MPNCIOFields
{

public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        // use default fields from the MPNC model
        Dumux::MPNCIOFields::initOutputModule(out);

        using VolumeVariables = typename OutputModule::VolumeVariables;
        using FS = typename VolumeVariables::FluidSystem;

        CommonIOFields::initOutputModule(out);

        out.addVolumeVariable([](const auto& v)
        { return v.capillaryPressure(); }, "pc");

        out.addVolumeVariable([](const auto& v){ return v.priVars().state(); },
                              "phasePriVarIdx");

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
};

} // end namespace Dumux::PoreNetwork

#endif
