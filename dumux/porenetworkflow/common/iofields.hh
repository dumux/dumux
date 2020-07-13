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
 * \ingroup PoreNetworkFlow
 * \copydoc Dumux::PNMCommonIOFields
 */
#ifndef DUMUX_PNM_COMMON_OUTPUT_FIELDS_HH
#define DUMUX_PNM_COMMON_OUTPUT_FIELDS_HH


namespace Dumux
{

/*!
 * \ingroup PoreNetworkFlow
 * \brief Adds output fields specific to all pore-network models
 */
class PNMCommonIOFields
{

public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        out.addField(out.problem().gridGeometry().coordinationNumber(), "coordinationNumber", OutputModule::FieldType::vertex);

        out.addField(out.problem().gridGeometry().poreLabel(), "poreLabel", OutputModule::FieldType::vertex);

        out.addVolumeVariable([](const auto& volVars){ return volVars.poreRadius(); }, "poreRadius");

        out.addField(out.problem().gridGeometry().throatLabel(), "throatLabel", OutputModule::FieldType::element);

        out.addFluxVariable([](const auto& fluxVars, const auto& fluxVarsCache)
                             { return fluxVarsCache.throatRadius(); }, "throatRadius");

        out.addFluxVariable([](const auto& fluxVars, const auto& fluxVarsCache)
                             { return fluxVarsCache.throatLength(); }, "throatLength");
    }
};

} // end namespace Dumux

#endif
