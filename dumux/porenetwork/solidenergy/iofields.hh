// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PNMSolidEnergyModel
 * \copydoc Dumux::PoreNetwork::SolidEnergyIOFields
 */
#ifndef DUMUX_PNM_SOLID_ENERGY_IO_FIELDS_HH
#define DUMUX_PNM_SOLID_ENERGY_IO_FIELDS_HH

#include <dumux/porenetwork/common/iofields.hh>
#include <dumux/io/name.hh>

namespace Dumux::PoreNetwork {

/*!
 * \ingroup PNMSolidEnergyModel
 * \brief Adds output fields specific to the PNM solid-energy model
 */
class SolidEnergyIOFields
{
public:

    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        CommonIOFields::initOutputModule(out);
        out.addVolumeVariable( [](const auto& v){ return v.temperature(); }, IOName::temperature());
    }

    template <class ModelTraits, class FluidSystem = void, class SolidSystem = void>
    static std::string primaryVariableName(int pvIdx, int state = 0)
    {
        return IOName::temperature();
    }
};

} // end namespace Dumux::PoreNetwork

#endif
