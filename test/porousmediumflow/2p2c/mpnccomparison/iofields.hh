// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPTwoCTests
 * \brief Adds I/O fields specific to the twop model
 */
#ifndef DUMUX_TWOPTWOC_MPNC_IO_FIELDS_HH
#define DUMUX_TWOPTWOC_MPNC_IO_FIELDS_HH

#include <dumux/io/name.hh>

namespace Dumux {

/*!
 * \ingroup TwoPTwoCTests
 * \brief Adds I/O fields specific to the two-phase two-component model
 */
class TwoPTwoCMPNCIOFields
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        using VolumeVariables = typename OutputModule::VolumeVariables;
        using FS = typename VolumeVariables::FluidSystem;

        // register standardized output fields
        out.addVolumeVariable([](const auto& v){ return v.porosity(); },
                              IOName::porosity());

        out.addVolumeVariable([](const auto& v){ return v.saturation(FS::phase0Idx); },
                              IOName::saturation<FS>(FS::phase0Idx));
        out.addVolumeVariable([](const auto& v){ return v.saturation(FS::phase1Idx); },
                              IOName::saturation<FS>(FS::phase1Idx));

        out.addVolumeVariable([](const auto& v){ return v.pressure(FS::phase0Idx); },
                              IOName::pressure<FS>(FS::phase0Idx));
        out.addVolumeVariable([](const auto& v){ return v.pressure(FS::phase1Idx); },
                              IOName::pressure<FS>(FS::phase1Idx));

        out.addVolumeVariable([](const auto& v){ return v.density(FS::phase0Idx); },
                              IOName::density<FS>(FS::phase0Idx));
        out.addVolumeVariable([](const auto& v){ return v.density(FS::phase1Idx); },
                              IOName::density<FS>(FS::phase1Idx));

        out.addVolumeVariable([](const auto& v){ return v.mobility(FS::phase0Idx); },
                              IOName::mobility<FS>(FS::phase0Idx));
        out.addVolumeVariable([](const auto& v){ return v.mobility(FS::phase1Idx); },
                              IOName::mobility<FS>(FS::phase1Idx));

        for (int i = 0; i < VolumeVariables::numFluidPhases(); ++i)
            for (int j = 0; j < VolumeVariables::numFluidComponents(); ++j)
                out.addVolumeVariable([i,j](const auto& v){ return v.massFraction(i,j); },
                                      IOName::massFraction<FS>(i, j));

        for (int i = 0; i < VolumeVariables::numFluidPhases(); ++i)
            for (int j = 0; j < VolumeVariables::numFluidComponents(); ++j)
                out.addVolumeVariable([i,j](const auto& v){ return v.moleFraction(i,j); },
                                      IOName::moleFraction<FS>(i, j));
    }
};

} // end namespace Dumux

#endif
