// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FluidStates
 * \brief Adapter class for fluid states with different indices
 */
#ifndef DUMUX_MATERIAL_FLUID_STATE_ADAPTER_HH
#define DUMUX_MATERIAL_FLUID_STATE_ADAPTER_HH

#include <dune/common/exceptions.hh>

namespace Dumux {

/*!
 * \ingroup FluidStates
 * \tparam FluidState the fluid state that is adapted to another fluid state's index set
 * \tparam AdapterPolicy
 * \brief Adapter class for fluid states with different indices
 * \note This is useful when writing a fluid system that has another fluid system as an ingredient
 *       To forward to the other fluid system (which variables to be sorted according to its own indices)
 *       create the fluid state adapter of the fluid state to forward and pass the adapter to this fluid state.
 *       You need to provide an adapter policy which maps the phase and component indices from the embedded fluid
 *       system's indices to the actual fluid system's indices.
 */
template<class FluidState, class AdapterPolicy>
class FluidStateAdapter
{
    using FluidSystem = typename AdapterPolicy::FluidSystem;
public:
    //! export scalar type
    using Scalar = typename FluidState::Scalar;

    //! export number of phases and components of the embedded fluid system
    static constexpr int numPhases = FluidSystem::numPhases;
    static constexpr int numComponents = FluidSystem::numComponents;

    /*!
     * \ingroup FluidStates
     * \brief Adapter class for fluid states with different indices
     * \param fluidState the fluid state to be wrapped
     */
    FluidStateAdapter(const FluidState& fluidState)
    : fluidState_(fluidState)
    {}

    /*!
     * \name Member functions forwarding to the fluid state mapping the indices
     */
    // \{
    int wettingPhase() const
    { DUNE_THROW(Dune::InvalidStateException, "wetting phase index cannot be mapped"); }

    Scalar moleFraction(int phaseIdx, int compIdx) const
    {
        Scalar sumMoleFrac = 0.0;
        for (int cIdx = 0; cIdx < numComponents; ++cIdx)
            sumMoleFrac += fluidState_.moleFraction(AdapterPolicy::phaseIdx(phaseIdx), AdapterPolicy::compIdx(cIdx));

        return fluidState_.moleFraction(AdapterPolicy::phaseIdx(phaseIdx), AdapterPolicy::compIdx(compIdx))/sumMoleFrac;
    }

    Scalar massFraction(int phaseIdx, int compIdx) const
    {
        Scalar sumMassFrac = 0.0;
        for (int cIdx = 0; cIdx < numComponents; ++cIdx)
            sumMassFrac += fluidState_.massFraction(AdapterPolicy::phaseIdx(phaseIdx), AdapterPolicy::compIdx(cIdx));

        return fluidState_.massFraction(AdapterPolicy::phaseIdx(phaseIdx), AdapterPolicy::compIdx(compIdx))/sumMassFrac;
    }

    Scalar averageMolarMass(int phaseIdx) const
    {
        Scalar sumXi = 0.0;
        Scalar sumXiMi = 0.0;
        for (int cIdx = 0; cIdx < numComponents; ++cIdx)
        {
            const auto xi = fluidState_.moleFraction(AdapterPolicy::phaseIdx(phaseIdx), AdapterPolicy::compIdx(cIdx));
            sumXi += xi;
            sumXiMi += xi*FluidSystem::molarMass(cIdx);
        }

        return sumXiMi/sumXi;
    }

    Scalar pressure(int phaseIdx) const
    { return fluidState_.pressure(AdapterPolicy::phaseIdx(phaseIdx)); }

    Scalar partialPressure(int phaseIdx, int compIdx) const
    {
        assert(FluidSystem::isGas(phaseIdx));
        return pressure(phaseIdx)*moleFraction(phaseIdx, compIdx);
    }

    Scalar temperature(int phaseIdx) const
    { return fluidState_.temperature(AdapterPolicy::phaseIdx(phaseIdx)); }

    Scalar temperature() const
    { return fluidState_.temperature(0); }
    // \}

private:
    const FluidState& fluidState_;
};

} // end namespace Dumux

#endif
