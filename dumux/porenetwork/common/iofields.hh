// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkModels
 * \copydoc Dumux::PoreNetwork::CommonIOFields
 */
#ifndef DUMUX_PNM_COMMON_OUTPUT_FIELDS_HH
#define DUMUX_PNM_COMMON_OUTPUT_FIELDS_HH

#include <dumux/io/vtk/function.hh>
#include <dumux/io/vtk/fieldtype.hh>

namespace Dumux::PoreNetwork {

/*!
 * \ingroup PoreNetworkModels
 * \brief Adds output fields specific to all pore-network models
 */
class CommonIOFields
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        out.addField(out.problem().gridGeometry().coordinationNumber(), "coordinationNumber", Vtk::FieldType::vertex);

        out.addField(out.problem().gridGeometry().poreLabel(), "poreLabel", Vtk::FieldType::vertex);

        out.addVolumeVariable([](const auto& volVars){ return volVars.poreInscribedRadius(); }, "poreInscribedRadius");

        out.addField(out.problem().gridGeometry().throatLabel(), "throatLabel", Vtk::FieldType::element);

        out.addFluxVariable([](const auto& fluxVars, const auto& fluxVarsCache)
                             { return fluxVarsCache.throatInscribedRadius(); }, "throatInscribedRadius");

        out.addFluxVariable([](const auto& fluxVars, const auto& fluxVarsCache)
                             { return fluxVarsCache.throatLength(); }, "throatLength");
    }
};

} // end namespace Dumux::PoreNetwork

#endif
