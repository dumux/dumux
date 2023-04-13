// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPNCModel
 * \brief Adds I/O fields specific to the twop-nc model.
 */

#ifndef DUMUX_TWOP_NC_IO_FIELDS_HH
#define DUMUX_TWOP_NC_IO_FIELDS_HH

#include <dumux/porousmediumflow/2p/iofields.hh>
#include <dumux/io/name.hh>

namespace Dumux
{

/*!
 * \ingroup TwoPNCModel
 * \brief Adds I/O fields specific to the TwoPNC model.
 */
class TwoPNCIOFields
{
public:
    template <class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        using VolumeVariables = typename OutputModule::VolumeVariables;
        using FluidSystem = typename VolumeVariables::FluidSystem;

        // use default fields from the 2p model
        TwoPIOFields::initOutputModule(out);

        // output additional to TwoP output:
        for (int phaseIdx = 0; phaseIdx < VolumeVariables::numFluidPhases(); ++phaseIdx)
        {
            for (int compIdx = 0; compIdx < VolumeVariables::numFluidComponents(); ++compIdx)
            {
                out.addVolumeVariable([phaseIdx,compIdx](const auto& v){ return v.moleFraction(phaseIdx,compIdx); },
                                      IOName::moleFraction<FluidSystem>(phaseIdx, compIdx));
                if (VolumeVariables::numFluidComponents() < 3)
                    out.addVolumeVariable([phaseIdx,compIdx](const auto& v){ return v.massFraction(phaseIdx,compIdx); },
                                          IOName::massFraction<FluidSystem>(phaseIdx, compIdx));
            }

            out.addVolumeVariable([phaseIdx](const auto& v){ return v.molarDensity(phaseIdx); },
                                    IOName::molarDensity<FluidSystem>(phaseIdx));
        }

        out.addVolumeVariable([](const auto& v){ return v.priVars().state(); },
                              IOName::phasePresence());
    }

    template <class ModelTraits, class FluidSystem, class SolidSystem = void>
    static std::string primaryVariableName(int pvIdx, int state)
    {
        using Indices = typename ModelTraits::Indices;
        static constexpr auto numStates = 3;
        using StringVec = std::array<std::string, numStates>;

        int idxSecComps;
        if (state == Indices::firstPhaseOnly
            || (state == Indices::bothPhases && ModelTraits::setMoleFractionsForFirstPhase()))
            idxSecComps = FluidSystem::phase0Idx;
        else
            idxSecComps = FluidSystem::phase1Idx;

        if (pvIdx > 1)
            return ModelTraits::useMoles() ? IOName::moleFraction<FluidSystem>(idxSecComps, pvIdx)
                                           : IOName::massFraction<FluidSystem>(idxSecComps, pvIdx);

        static const StringVec p0s1SwitchedPvNames = {
            ModelTraits::useMoles() ? IOName::moleFraction<FluidSystem>(FluidSystem::phase0Idx, FluidSystem::comp1Idx)
                                    : IOName::massFraction<FluidSystem>(FluidSystem::phase0Idx, FluidSystem::comp1Idx),
            ModelTraits::useMoles() ? IOName::moleFraction<FluidSystem>(FluidSystem::phase1Idx, FluidSystem::comp0Idx)
                                    : IOName::massFraction<FluidSystem>(FluidSystem::phase1Idx, FluidSystem::comp0Idx),
            IOName::saturation<FluidSystem>(FluidSystem::phase1Idx)};

        static const StringVec p1s0SwitchedPvNames = {
            ModelTraits::useMoles() ? IOName::moleFraction<FluidSystem>(FluidSystem::phase0Idx, FluidSystem::comp1Idx)
                                    : IOName::massFraction<FluidSystem>(FluidSystem::phase0Idx, FluidSystem::comp1Idx),
            ModelTraits::useMoles() ? IOName::moleFraction<FluidSystem>(FluidSystem::phase1Idx, FluidSystem::comp0Idx)
                                    : IOName::massFraction<FluidSystem>(FluidSystem::phase1Idx, FluidSystem::comp0Idx),
            IOName::saturation<FluidSystem>(FluidSystem::phase0Idx)};

        switch (ModelTraits::priVarFormulation())
        {
        case TwoPFormulation::p0s1:
            return pvIdx == 0 ? IOName::pressure<FluidSystem>(FluidSystem::phase0Idx)
                              : p0s1SwitchedPvNames[state-1];
        case TwoPFormulation::p1s0:
            return pvIdx == 0 ? IOName::pressure<FluidSystem>(FluidSystem::phase1Idx)
                              : p1s0SwitchedPvNames[state-1];
        default: DUNE_THROW(Dune::InvalidStateException, "Invalid formulation ");
        }
    }
};


} // end namespace Dumux

#endif
