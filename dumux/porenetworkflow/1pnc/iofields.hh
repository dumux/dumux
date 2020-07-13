// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup PoreNetworkOnePNCModel
 * \copydoc Dumux::PNMOnePNCIOFields
 */
#ifndef DUMUX_PNM_ONEP_NC_IO_FIELDS_HH
#define DUMUX_PNM_ONEP_NC_IO_FIELDS_HH

#include <dumux/porousmediumflow/1pnc/iofields.hh>
#include <dumux/porenetworkflow/common/iofields.hh>

namespace Dumux
{

/*!
 * \ingroup PoreNetworkOnePNCModel
 * \brief Adds output fields specific to the PNM 1pnc model
 */
template<class FluidSystem>
class PNMOnePNCIOFields
{

public:
    template<class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        OnePNCIOFields::initOutputModule(out);
        PNMCommonIOFields::initOutputModule(out);

        out.addFluxVariable([](const auto& fluxVars, const auto& fluxVarsCache)
                              { return fluxVarsCache.transmissibility(0); }, "transmissibility");

        auto volumeFlux = [](const auto& fluxVars, const auto& fluxVarsCache)
        {
            auto upwindTerm = [](const auto& volVars) { return volVars.mobility(0); };
            using std::abs;
            return abs(fluxVars.advectiveFlux(0, upwindTerm));
        };
        out.addFluxVariable(volumeFlux, "volumeFlux");
    }
};

} // end namespace Dumux

#endif // DUMUX_PNM_ONEP_NC_VTK_OUTPUT_FIELDS_HH
