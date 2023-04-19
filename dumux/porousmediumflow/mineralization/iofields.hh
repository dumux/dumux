// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MineralizationModel
 * \brief Adds I/O fields specific to the models considering
 *        mineralization processes.
 */

#ifndef DUMUX_MINERALIZATION_IO_FIELDS_HH
#define DUMUX_MINERALIZATION_IO_FIELDS_HH

#include <dumux/io/name.hh>

namespace Dumux {

/*!
 * \ingroup MineralizationModel
 * \brief Adds I/O fields specific to a NCMin model.
 */
template<class NonMineralizationIOFields>
class MineralizationIOFields
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        using SolidSystem = typename OutputModule::VolumeVariables::SolidSystem;

        // output of the model without mineralization
        NonMineralizationIOFields::initOutputModule(out);

        // additional output
        for (int i = 0; i < SolidSystem::numComponents - SolidSystem::numInertComponents; ++i)
        {
            out.addVolumeVariable([i](const auto& v){ return v.solidVolumeFraction(i); },
                                  IOName::solidVolumeFraction<SolidSystem>(i));
        }
    }

    template <class ModelTraits, class FluidSystem, class SolidSystem>
    static std::string primaryVariableName(int pvIdx, int state = 0)
    {
        static constexpr int nonMinNumEq = ModelTraits::numEq() - ModelTraits::numSolidComps() + ModelTraits::numInertSolidComps();

        if (pvIdx < nonMinNumEq)
            return NonMineralizationIOFields::template primaryVariableName<ModelTraits, FluidSystem, SolidSystem>(pvIdx, state);
        else
            return IOName::solidVolumeFraction<SolidSystem>(pvIdx - nonMinNumEq);
    }
};

} // end namespace Dumux

#endif
