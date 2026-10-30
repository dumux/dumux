// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
/*!
 * \file
 * \ingroup TwoPTwoCTests
 * \brief Water-air fluid system with the fluid properties of the heat-pipe solution of Huang et al. (2015).
 */
#ifndef DUMUX_HEATPIPE_REFERENCE_FLUID_SYSTEM_HH
#define DUMUX_HEATPIPE_REFERENCE_FLUID_SYSTEM_HH

#include <cmath>

#include <dumux/material/fluidsystems/h2oair.hh>

namespace Dumux::FluidSystems {

/*!
 * \ingroup TwoPTwoCTests
 * \brief Water-air fluid system with the fluid properties of the semi-analytical heat-pipe
 *        solution of Udell and Fitch (1985) as extended by Huang et al. (2015), see
 *        test_heatpipe_odesolver.cc. All values are taken from the supplementary code of Huang et al.
 *
 * - constant liquid density and viscosity, (almost) no air dissolved in the liquid
 * - ideal gas mixture with a mole-fraction weighted viscosity
 * - constant diffusion coefficient of water vapor in air
 * - Clausius-Clapeyron vapor pressure with a constant latent heat and the Kelvin correction
 * - a common sensible enthalpy for all components, plus the latent heat for water vapor
 */
template<class Scalar>
class HeatPipeReferenceFluidSystem
: public H2OAir<Scalar>
{
    using Parent = H2OAir<Scalar>;

public:
    using typename Parent::ParameterCache;
    using Parent::liquidPhaseIdx;
    using Parent::gasPhaseIdx;
    using Parent::H2OIdx;
    using Parent::AirIdx;

    static constexpr Scalar waterMolarMass = 0.018016; // kg/mol
    static constexpr Scalar airMolarMass = 0.028; // kg/mol
    static constexpr Scalar gasConstant = 8.314; // J/(mol K)
    static constexpr Scalar liquidDensity = 1000.0; // kg/m^3
    static constexpr Scalar liquidViscosity = 2.938e-4; // Pa s
    static constexpr Scalar airViscosity = 2.094e-5; // Pa s
    static constexpr Scalar vaporViscosity = 1.2e-5; // Pa s
    static constexpr Scalar gasDiffusionCoefficient = 2.6e-6; // m^2/s, macroscopic (porous medium) value
    static constexpr Scalar latentHeat = 2.258e6; // J/kg
    static constexpr Scalar boilingTemperature = 373.15; // K, at the boiling pressure
    static constexpr Scalar boilingPressure = 101325.0; // Pa
    static constexpr Scalar specificHeatCapacity = 4187.0; // J/(kg K), only affects the transient

    static constexpr Scalar molarMass(int compIdx)
    { return compIdx == H2OIdx ? waterMolarMass : airMolarMass; }

    /*!
     * \brief Vapor pressure from the Clausius-Clapeyron equation with the Kelvin correction
     *        \f$p_v = p_\mathrm{sat}(T)\exp[-(p_v - p_w)M^w/(\rho_w R T)]\f$.
     *
     * This is the form \f$p_c - x_g^a p_g = p_v - p_w\f$ of Huang et al. (2015).
     * Since the exponent is small, a fixed-point iteration converges within a few steps.
     */
    template<class FluidState>
    static Scalar vaporPressure(const FluidState& state, int compIdx)
    {
        using std::exp;
        const Scalar t = state.temperature(liquidPhaseIdx);
        const Scalar pw = state.pressure(liquidPhaseIdx);
        const Scalar pSat = boilingPressure*exp((1.0/boilingTemperature - 1.0/t)*latentHeat*waterMolarMass/gasConstant);
        Scalar pv = pSat;
        for (int i = 0; i < 10; ++i)
            pv = pSat*exp(-(pv - pw)*waterMolarMass/(liquidDensity*gasConstant*t));
        return pv;
    }

    template<class FluidState>
    static Scalar fugacityCoefficient(const FluidState& state, int phaseIdx, int compIdx)
    {
        if (phaseIdx == gasPhaseIdx)
            return 1.0;

        // A finite, large Henry constant makes dissolved air negligible while
        // keeping the compositional equilibrium solver well-defined.
        const Scalar p = state.pressure(liquidPhaseIdx);
        return compIdx == H2OIdx ? vaporPressure(state, H2OIdx)/p : 1e20/p;
    }

    template<class FluidState>
    static Scalar fugacityCoefficient(const FluidState& state, const ParameterCache&, int phaseIdx, int compIdx)
    { return fugacityCoefficient(state, phaseIdx, compIdx); }

    template<class FluidState>
    static Scalar molarDensity(const FluidState& state, int phaseIdx)
    {
        if (phaseIdx == liquidPhaseIdx)
            return liquidDensity/state.averageMolarMass(liquidPhaseIdx);
        return state.pressure(gasPhaseIdx)/(gasConstant*state.temperature(gasPhaseIdx));
    }

    template<class FluidState>
    static Scalar molarDensity(const FluidState& state, const ParameterCache&, int phaseIdx)
    { return molarDensity(state, phaseIdx); }

    template<class FluidState>
    static Scalar density(const FluidState& state, int phaseIdx)
    { return molarDensity(state, phaseIdx)*state.averageMolarMass(phaseIdx); }

    template<class FluidState>
    static Scalar density(const FluidState& state, const ParameterCache&, int phaseIdx)
    { return density(state, phaseIdx); }

    template<class FluidState>
    static Scalar viscosity(const FluidState& state, int phaseIdx)
    {
        if (phaseIdx == liquidPhaseIdx)
            return liquidViscosity;
        return state.moleFraction(gasPhaseIdx, AirIdx)*airViscosity
               + state.moleFraction(gasPhaseIdx, H2OIdx)*vaporViscosity;
    }

    template<class FluidState>
    static Scalar viscosity(const FluidState& state, const ParameterCache&, int phaseIdx)
    { return viscosity(state, phaseIdx); }

    template<class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState& state, int phaseIdx, int compIIdx, int compJIdx)
    {
        if (phaseIdx == gasPhaseIdx)
            return gasDiffusionCoefficient;
        return Parent::binaryDiffusionCoefficient(state, phaseIdx, compIIdx, compJIdx);
    }

    template<class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState& state, const ParameterCache&,
                                             int phaseIdx, int compIIdx, int compJIdx)
    { return binaryDiffusionCoefficient(state, phaseIdx, compIIdx, compJIdx); }

    template<class FluidState>
    static Scalar componentEnthalpy(const FluidState& state, int phaseIdx, int compIdx)
    {
        // A common sensible enthalpy cancels from the steady advective energy
        // flux when the total component mass flux is zero. Only water's fixed
        // latent heat remains, as assumed by the reference energy balance.
        const Scalar sensible = specificHeatCapacity*(state.temperature(phaseIdx) - 273.15);
        return sensible + ((phaseIdx == gasPhaseIdx && compIdx == H2OIdx) ? latentHeat : 0.0);
    }

    template<class FluidState>
    static Scalar enthalpy(const FluidState& state, int phaseIdx)
    {
        return componentEnthalpy(state, phaseIdx, H2OIdx)*state.massFraction(phaseIdx, H2OIdx)
             + componentEnthalpy(state, phaseIdx, AirIdx)*state.massFraction(phaseIdx, AirIdx);
    }

    template<class FluidState>
    static Scalar enthalpy(const FluidState& state, const ParameterCache&, int phaseIdx)
    { return enthalpy(state, phaseIdx); }

    template<class FluidState>
    static Scalar heatCapacity(const FluidState&, int)
    { return specificHeatCapacity; }

    template<class FluidState>
    static Scalar heatCapacity(const FluidState& state, const ParameterCache&, int phaseIdx)
    { return heatCapacity(state, phaseIdx); }
};

} // namespace Dumux::FluidSystems

#endif
