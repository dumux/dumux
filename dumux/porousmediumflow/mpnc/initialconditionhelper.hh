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
 * \ingroup MPNCModel
 * \brief A helper function to get the correct initial conditions by updating the fluidstate and defining the primary variables needed for equilibrium mpnc models for the MPNC model
 */
#ifndef DUMUX_MPNC_INITIALCONDITION_HELPER_HH
#define DUMUX_MPNC_INITIALCONDITION_HELPER_HH

#include <dumux/material/constraintsolvers/misciblemultiphasecomposition.hh>
#include <dumux/material/constraintsolvers/computefromreferencephase.hh>

namespace Dumux {

namespace MPNCInitialConditions{

    struct AllPhasesPresent { int refPhaseIdx; };
    struct NotAllPhasesPresent { int refPhaseIdx; };
}

template <class Scalar, class PrimaryVariables, class FluidSystem, class ModelTraits>
struct MPNCInitialConditionHelper
{

    template<class FluidState, class ParamCache>
    static auto solveForPrimaryVariables(FluidState& fluidState,
                                         ParamCache& paramCache,
                                         const MPNCInitialConditions::AllPhasesPresent& allPhases)
    {

        // make the fluid state consistent with local thermodynamic equilibrium
        using MiscibleMultiPhaseComposition = Dumux::MiscibleMultiPhaseComposition<Scalar, FluidSystem>;

        MiscibleMultiPhaseComposition::solve(fluidState, paramCache, allPhases.refPhaseIdx);

        static constexpr auto numComponents = FluidSystem::numComponents;
        static constexpr auto numPhases = FluidSystem::numPhases;
        static constexpr auto fug0Idx = ModelTraits::Indices::fug0Idx;
        static constexpr auto s0Idx = ModelTraits::Indices::s0Idx;
        static constexpr auto p0Idx = ModelTraits::Indices::p0Idx;

        PrimaryVariables values(0.0);
        // all N component fugacities
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            values[fug0Idx + compIdx] = fluidState.fugacity(0, compIdx);

        // first M - 1 saturations
        for (int phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx)
            values[s0Idx + phaseIdx] = fluidState.saturation(phaseIdx);

        static constexpr auto pressureFormulation = ModelTraits::pressureFormulation();

        // first pressure
        if(pressureFormulation == MpNcPressureFormulation::mostWettingFirst)
            values[p0Idx] = fluidState.pressure(/*phaseIdx=*/0);
        else if(pressureFormulation == MpNcPressureFormulation::leastWettingFirst)
            values[p0Idx] = fluidState.pressure(numPhases-1);
        else
            DUNE_THROW(Dune::InvalidStateException,"unknown pressure formulation");

        return values;
    }

    template<class FluidState, class ParamCache>
    static auto solveForPrimaryVariables(FluidState& fluidState,
                                         ParamCache& paramCache,
                                         const MPNCInitialConditions::NotAllPhasesPresent& notAllPhases)
    {

        // make the fluid state consistent with local thermodynamic equilibrium
        using ComputeFromReferencePhase = ComputeFromReferencePhase<Scalar, FluidSystem>;

        ComputeFromReferencePhase::solve(fluidState,
                                     paramCache,
                                     notAllPhases.refPhaseIdx);

        static constexpr auto numComponents = FluidSystem::numComponents;
        static constexpr auto numPhases = FluidSystem::numPhases;
        static constexpr auto fug0Idx = ModelTraits::Indices::fug0Idx;
        static constexpr auto s0Idx = ModelTraits::Indices::s0Idx;
        static constexpr auto p0Idx = ModelTraits::Indices::p0Idx;

        PrimaryVariables values(0.0);
        // all N component fugacities
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            values[fug0Idx + compIdx] = fluidState.fugacity(0, compIdx);

        // first M - 1 saturations
        for (int phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx)
            values[s0Idx + phaseIdx] = fluidState.saturation(phaseIdx);

        static constexpr auto pressureFormulation = ModelTraits::pressureFormulation();
        // first pressure
        if(pressureFormulation == MpNcPressureFormulation::mostWettingFirst)
            values[p0Idx] = fluidState.pressure(/*phaseIdx=*/0);
        else if(pressureFormulation == MpNcPressureFormulation::leastWettingFirst)
            values[p0Idx] = fluidState.pressure(numPhases-1);
        else
            DUNE_THROW(Dune::InvalidStateException,"unknown pressure formulation");

        return values;
    }
};


} // end namespace Dumux

#endif
