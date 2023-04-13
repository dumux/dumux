// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MPNCModel
 * \brief A helper function to get the correct initial conditions by updating the fluidstate and defining the primary variables needed for equilibrium mpnc models for the MPNC model
 */
#ifndef DUMUX_MPNC_INITIALCONDITION_HELPER_HH
#define DUMUX_MPNC_INITIALCONDITION_HELPER_HH

#include <dumux/material/constraintsolvers/misciblemultiphasecomposition.hh>
#include <dumux/material/constraintsolvers/computefromreferencephase.hh>

#include "pressureformulation.hh"

namespace Dumux {

namespace MPNCInitialConditions {

struct AllPhasesPresent { int refPhaseIdx; };
struct NotAllPhasesPresent { int refPhaseIdx; };

} // namespace MPNCInitialConditions

template<class PrimaryVariables, class FluidSystem, class ModelTraits>
class MPNCInitialConditionHelper
{
    using Scalar = typename PrimaryVariables::value_type;
    using AllPhasesPresent = MPNCInitialConditions::AllPhasesPresent;
    using NotAllPhasesPresent = MPNCInitialConditions::NotAllPhasesPresent;

public:
    template<class FluidState>
    void solve(FluidState& fs, const AllPhasesPresent& allPhases) const
    {
        typename FluidSystem::ParameterCache paramCache;
        MiscibleMultiPhaseComposition<Scalar, FluidSystem>::solve(fs,
                                                                  paramCache,
                                                                  allPhases.refPhaseIdx);
    }

    template<class FluidState>
    void solve(FluidState& fs, const NotAllPhasesPresent& notAllPhases) const
    {
        typename FluidSystem::ParameterCache paramCache;
        ComputeFromReferencePhase<Scalar, FluidSystem>::solve(fs,
                                                              paramCache,
                                                              notAllPhases.refPhaseIdx);
    }

    template<class FluidState>
    PrimaryVariables getPrimaryVariables(const FluidState& fs) const
    {
        PrimaryVariables priVars(0.0);

        static constexpr auto numComponents = FluidSystem::numComponents;
        static constexpr auto numPhases = FluidSystem::numPhases;
        static constexpr auto fug0Idx = ModelTraits::Indices::fug0Idx;
        static constexpr auto s0Idx = ModelTraits::Indices::s0Idx;
        static constexpr auto p0Idx = ModelTraits::Indices::p0Idx;

        // all N component fugacities
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            priVars[fug0Idx + compIdx] = fs.fugacity(0, compIdx);

        // first M - 1 saturations
        for (int phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx)
            priVars[s0Idx + phaseIdx] = fs.saturation(phaseIdx);

        static constexpr auto pressureFormulation = ModelTraits::pressureFormulation();
        if (pressureFormulation == MpNcPressureFormulation::mostWettingFirst)
            priVars[p0Idx] = fs.pressure(/*phaseIdx=*/0);
        else if (pressureFormulation == MpNcPressureFormulation::leastWettingFirst)
            priVars[p0Idx] = fs.pressure(numPhases-1);
        else
            DUNE_THROW(Dune::InvalidStateException,"unknown pressure formulation");

        return priVars;
    }
};

} // end namespace Dumux

#endif
