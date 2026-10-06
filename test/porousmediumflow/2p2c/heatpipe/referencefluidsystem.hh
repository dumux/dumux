// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef DUMUX_HEATPIPE_REFERENCE_FLUID_SYSTEM_HH
#define DUMUX_HEATPIPE_REFERENCE_FLUID_SYSTEM_HH

#include <dumux/material/fluidsystems/h2oair.hh>

namespace Dumux::FluidSystems {

// Match the density assumptions of test_heatpipe_odesolver.cc without changing
// the local component viscosities or Wilke's gas-mixture viscosity rule.
struct HeatPipeReferencePolicy : H2OAirDefaultPolicy<>
{
    static constexpr bool useH2ODensityAsLiquidMixtureDensity() { return true; }
    static constexpr bool useIdealGasDensity() { return true; }
};

template<class Scalar>
class HeatPipeReferenceFluidSystem
: public H2OAir<Scalar, Components::TabulatedComponent<Components::H2O<Scalar>>, HeatPipeReferencePolicy>
{
    using Parent = H2OAir<Scalar, Components::TabulatedComponent<Components::H2O<Scalar>>, HeatPipeReferencePolicy>;
public:
    using typename Parent::ParameterCache;

    template<class FluidState>
    static Scalar fugacityCoefficient(const FluidState& state, int phaseIdx, int compIdx)
    {
        // A finite, large Henry constant makes dissolved air negligible while
        // keeping the compositional equilibrium solver well-defined.
        if (phaseIdx == Parent::liquidPhaseIdx && compIdx == Parent::AirIdx)
            return Scalar(1e20)/state.pressure(phaseIdx);
        return Parent::fugacityCoefficient(state, phaseIdx, compIdx);
    }

    template<class FluidState>
    static Scalar fugacityCoefficient(const FluidState& state, const ParameterCache&, int phaseIdx, int compIdx)
    { return fugacityCoefficient(state, phaseIdx, compIdx); }

    template<class FluidState>
    static Scalar componentEnthalpy(const FluidState& state, int phaseIdx, int compIdx)
    {
        // A common sensible enthalpy cancels from the steady advective energy
        // flux when the total component mass flux is zero. Only water's fixed
        // latent heat remains, as assumed by the reference energy balance.
        // Retain positive heat capacity for the transient approach to equilibrium.
        const Scalar sensible = Scalar(4187)*(state.temperature(phaseIdx) - Scalar(273.15));
        return sensible + ((phaseIdx == Parent::gasPhaseIdx && compIdx == Parent::H2OIdx)
                           ? Scalar(2.258e6) : Scalar(0));
    }

    template<class FluidState>
    static Scalar enthalpy(const FluidState& state, int phaseIdx)
    {
        return componentEnthalpy(state, phaseIdx, Parent::H2OIdx)*state.massFraction(phaseIdx, Parent::H2OIdx)
             + componentEnthalpy(state, phaseIdx, Parent::AirIdx)*state.massFraction(phaseIdx, Parent::AirIdx);
    }

    template<class FluidState>
    static Scalar enthalpy(const FluidState& state, const ParameterCache&, int phaseIdx)
    { return enthalpy(state, phaseIdx); }

    template<class FluidState>
    static Scalar heatCapacity(const FluidState&, int)
    { return Scalar(4187); }

    template<class FluidState>
    static Scalar heatCapacity(const FluidState& state, const ParameterCache&, int phaseIdx)
    { return heatCapacity(state, phaseIdx); }
};

} // namespace Dumux::FluidSystems

#endif
