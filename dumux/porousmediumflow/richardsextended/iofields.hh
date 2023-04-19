// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ExtendedRichardsModel
 * \brief Adds I/O fields specific to the extended Richards model.
 */

#ifndef DUMUX_RICHARDSEXTENDED_IO_FIELDS_HH
#define DUMUX_RICHARDSEXTENDED_IO_FIELDS_HH

#include <dumux/io/name.hh>
#include <dumux/porousmediumflow/richards/iofields.hh>

namespace Dumux {

/*!
 * \ingroup ExtendedRichardsModel
 * \brief Adds I/O fields specific to the extended Richards model.
 */
class ExtendedRichardsIOFields : public RichardsIOFields
{
    using ParentType = RichardsIOFields;
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        using VV = typename OutputModule::VolumeVariables;
        using FS = typename VV::FluidSystem;

        ParentType::initOutputModule(out);

        out.addVolumeVariable([](const auto& v){ return v.moleFraction(FS::phase1Idx, FS::comp0Idx); },
                                  IOName::moleFraction<FS>(FS::phase1Idx, FS::comp0Idx));
        out.addVolumeVariable([](const auto& v){ return v.priVars().state(); },
                                  IOName::phasePresence());
    }
};

    template<class ModelTraits, class FluidSystem, class SolidSystem = void>
    static std::string primaryVariableName(int pvIdx, int state)
    {
        using Indices = typename ModelTraits::Indices;

        if (state == Indices::gasPhaseOnly)
            return IOName::moleFraction<FluidSystem>(FluidSystem::phase1Idx, FluidSystem::phase0Idx);
        else
            return IOName::pressure<FluidSystem>(FluidSystem::phase0Idx);
    }
} // end namespace Dumux

#endif
