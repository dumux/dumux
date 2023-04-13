// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPOneCModel
 * \copydoc Dumux::TwoPOneCIOFields
 */

#ifndef DUMUX_TWOP_ONEC_IO_FIELDS_HH
#define DUMUX_TWOP_ONEC_IO_FIELDS_HH

#include <dumux/porousmediumflow/2p/iofields.hh>
#include <dumux/io/name.hh>

namespace Dumux {

/*!
 * \ingroup TwoPOneCModel
 * \brief Adds I/O fields specific to two-phase one-component model.
 */
class TwoPOneCIOFields
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        // use default fields from the 2p model
        TwoPIOFields::initOutputModule(out);

        // output additional to TwoP output:
        out.addVolumeVariable([](const auto& v){ return v.priVars().state(); },
                              IOName::phasePresence());
    }

    template <class ModelTraits, class FluidSystem, class SolidSystem = void>
    static std::string primaryVariableName(int pvIdx, int state)
    {
        using Indices = typename ModelTraits::Indices;

        if (ModelTraits::priVarFormulation() == TwoPFormulation::p0s1)
            return (pvIdx == 0) ? IOName::pressure<FluidSystem>(FluidSystem::phase0Idx) :
                                  (state == Indices::twoPhases)
                                  ? IOName::saturation<FluidSystem>(FluidSystem::phase1Idx)
                                  : IOName::temperature();
        else
            return (pvIdx == 0) ? IOName::pressure<FluidSystem>(FluidSystem::phase1Idx) :
                                  (state == Indices::twoPhases)
                                  ? IOName::saturation<FluidSystem>(FluidSystem::phase0Idx)
                                  : IOName::temperature();
    }
};

} // end namespace Dumux

#endif
