// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPTwoCTests
 * \brief Semi-analytical reference solution for the heatpipe benchmark.
 *
 * Integrates the steady-state ODE system of Udell and Fitch (1985), in the
 * formulation of Huang et al. (2015), from the left (Dirichlet) boundary
 * towards the heat source until the wetting phase dries out. Fluid properties,
 * capillary pressure and relative permeabilities are evaluated with the same
 * fluid system and constitutive laws as the numerical model.
 */
#include <config.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/math.hh>
#include <dumux/experimental/ode/odesolver.hh>
#include <dumux/material/constants.hh>
#include <dumux/material/fluidstates/compositional.hh>
#include <dumux/material/fluidsystems/h2oair.hh>
#include <dumux/material/fluidmatrixinteractions/2p/heatpipelaw.hh>

namespace Dumux {

/*!
 * \ingroup TwoPTwoCTests
 * \brief The heatpipe ODE system for the effective wetting-phase saturation,
 *        the gas-phase pressure, the gas-phase air mole fraction and the temperature.
 */
class HeatPipeReferenceODE
{
public:
    using Scalar = double;
    using SolutionVector = Dune::FieldVector<Scalar, 4>;
    using ResidualType = SolutionVector;
    using JacobianMatrix = Dune::FieldMatrix<Scalar, 4, 4>;
    using Variables = Experimental::ODEVariables<SolutionVector>;
    // same fluid system as the numerical model (see properties.hh)
    using FluidSystem = FluidSystems::H2OAir<Scalar>;

    enum Component
    {
        effectiveSaturationIdx,
        gasPressureIdx,
        gasAirMoleFractionIdx,
        temperatureIdx
    };

    HeatPipeReferenceODE()
    : permeability_(getParam<Scalar>("Problem.Permeability"))
    , heatFlux_(getParam<Scalar>("Problem.HeatFlux"))
    , lambdaSolid_(getParam<Scalar>("Component.SolidThermalConductivity"))
    , pcKrSwCurve_(PcKrSwCurve::Params(surfaceTension_, std::sqrt(porosity_/permeability_)),
                   PcKrSwCurve::EffToAbsParams(swr_, 0.0))
    {}

    SolutionVector initialState() const
    {
        SolutionVector result;
        result[effectiveSaturationIdx] = (swBc_ - swr_)/(1.0 - swr_);
        result[gasPressureIdx] = pgBc_;
        // local equilibrium, neglecting the air dissolved in the liquid phase
        result[gasAirMoleFractionIdx] = 1.0 - FluidSystem::H2O::vaporPressure(tBc_)/pgBc_;
        result[temperatureIdx] = tBc_;
        return result;
    }

    void rhs(const Variables& vars, ResidualType& rhs) const
    {
        const auto& z = vars.dofs();
        const auto se = z[effectiveSaturationIdx];
        const auto pg = z[gasPressureIdx];
        const auto xa = z[gasAirMoleFractionIdx];
        const auto temperature = z[temperatureIdx];

        const auto sw = saturation(se);
        const auto pc = pcKrSwCurve_.pc(sw);
        const auto dpcDse = capillaryPressureDerivative_(sw);
        const auto krg = pcKrSwCurve_.krn(sw);
        const auto krl = pcKrSwCurve_.krw(sw);

        const auto fluidState = fluidState_(pg, pg - pc, xa, temperature);
        const auto rhoW = FluidSystem::density(fluidState, liquidPhaseIdx);
        const auto muW = FluidSystem::viscosity(fluidState, liquidPhaseIdx);
        const auto muG = FluidSystem::viscosity(fluidState, gasPhaseIdx);
        const auto dPm = millingtonQuirk_(1.0 - sw, FluidSystem::binaryDiffusionCoefficient(fluidState, gasPhaseIdx, H2OIdx, AirIdx));
        const auto lambda = somerton_(sw, FluidSystem::thermalConductivity(fluidState, liquidPhaseIdx),
                                          FluidSystem::thermalConductivity(fluidState, gasPhaseIdx));

        // the ODE system is derived for an ideal gas phase
        const auto rhoG = (molarMassAir_*xa + molarMassWater_*(1.0 - xa))*pg/(gasConstant_*temperature);
        const auto nuG = muG/rhoG;
        const auto nuW = muW/rhoW;
        const auto beta = nuW/nuG;
        const auto alpha = 1.0 + pc/(rhoW*latentHeat_);
        const auto xi = (1.0/krg)*(1.0 + rhoW*gasConstant_*temperature
                                         /(pg*molarMassWater_*(1.0 - xa)))
                        + beta/krl;
        const auto delta = rhoW*latentHeat_*latentHeat_*permeability_*alpha/(lambda*nuG*temperature);
        const auto zeta = permeability_*rhoW*gasConstant_*temperature/(molarMassWater_*rhoG*nuG*dPm)
                          *xa/(1.0 - xa)
                          *(pg*molarMassWater_/(rhoW*gasConstant_*temperature) + 1.0/(1.0 - xa));
        const auto eta = delta/(delta + xi + zeta);

        rhs[effectiveSaturationIdx] = -(1.0/(1.0 - xa) + beta*krg/krl)
                                      *eta*heatFlux_*nuG/(permeability_*latentHeat_*krg*dpcDse);
        rhs[gasPressureIdx] = -(eta*heatFlux_*nuG/(permeability_*latentHeat_*krg))/(1.0 - xa);
        rhs[gasAirMoleFractionIdx] = eta*heatFlux_*xa/(latentHeat_*dPm*rhoG*(1.0 - xa));
        rhs[temperatureIdx] = -heatFlux_*(1.0 - eta)/lambda;
    }

    Scalar saturation(const Scalar se) const
    { return swr_ + se*(1.0 - swr_); }

    //! The temperature gradient in the dry zone, where heat is transported by conduction only
    Scalar dryTemperatureGradient() const
    { return -heatFlux_/somerton_(0.0, 0.0, airThermalConductivity_()); }

private:
    static constexpr int liquidPhaseIdx = FluidSystem::liquidPhaseIdx;
    static constexpr int gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr int H2OIdx = FluidSystem::H2OIdx;
    static constexpr int AirIdx = FluidSystem::AirIdx;
    using PcKrSwCurve = FluidMatrix::HeatPipeLaw<Scalar>;
    using FluidState = CompositionalFluidState<Scalar, FluidSystem>;

    FluidState fluidState_(const Scalar pg, const Scalar pw, const Scalar xa, const Scalar temperature) const
    {
        FluidState fluidState;
        fluidState.setTemperature(temperature);
        fluidState.setPressure(gasPhaseIdx, pg);
        fluidState.setPressure(liquidPhaseIdx, pw);
        fluidState.setMoleFraction(gasPhaseIdx, AirIdx, xa);
        fluidState.setMoleFraction(gasPhaseIdx, H2OIdx, 1.0 - xa);
        fluidState.setMoleFraction(liquidPhaseIdx, AirIdx, 0.0);
        fluidState.setMoleFraction(liquidPhaseIdx, H2OIdx, 1.0);
        return fluidState;
    }

    Scalar capillaryPressureDerivative_(const Scalar sw) const
    {
        // HeatPipeLaw::dpc_dsw can currently not be instantiated, so use a central difference
        const auto eps = 1e-8;
        const auto dpcDsw = (pcKrSwCurve_.pc(sw + eps) - pcKrSwCurve_.pc(sw - eps))/(2.0*eps);
        return dpcDsw*(1.0 - swr_);
    }

    //! ThermalConductivitySomertonTwoP (default of the TwoPTwoCNI model)
    Scalar somerton_(const Scalar sw, const Scalar lambdaLiquid, const Scalar lambdaGas) const
    {
        using std::pow;
        using std::sqrt;
        using std::max;
        const auto lambdaSaturated = lambdaSolid_*pow(lambdaLiquid/lambdaSolid_, porosity_);
        const auto lambdaDry = lambdaSolid_*pow(lambdaGas/lambdaSolid_, porosity_);
        return lambdaDry + sqrt(max(sw, 0.0))*(lambdaSaturated - lambdaDry);
    }

    //! DiffusivityMillingtonQuirk (default of the TwoPTwoC model)
    Scalar millingtonQuirk_(const Scalar sg, const Scalar binaryDiffusionCoefficient) const
    {
        using std::pow;
        using std::cbrt;
        using std::max;
        const auto positiveSg = max(sg, 0.0);
        return porosity_*pow(positiveSg, 3)*cbrt(porosity_*positiveSg)*binaryDiffusionCoefficient;
    }

    Scalar airThermalConductivity_() const
    { return FluidSystem::thermalConductivity(fluidState_(pgBc_, pgBc_, 1.0, tBc_), gasPhaseIdx); }

    // spatial parameters, see spatialparams.hh
    static constexpr Scalar porosity_ = 0.4;
    static constexpr Scalar swr_ = 0.15;
    static constexpr Scalar surfaceTension_ = 0.0588;
    // Dirichlet values at the left boundary, see problem.hh
    static constexpr Scalar pgBc_ = 1.013e5;
    static constexpr Scalar swBc_ = 0.99;
    static constexpr Scalar tBc_ = 341.75;
    // latent heat of vaporization at the normal boiling point, a parameter of the ODE system
    static constexpr Scalar latentHeat_ = 2.258e6;
    static constexpr Scalar gasConstant_ = Constants<Scalar>::R;
    static constexpr Scalar molarMassWater_ = FluidSystem::H2O::molarMass();
    static constexpr Scalar molarMassAir_ = FluidSystem::Air::molarMass();

    Scalar permeability_;
    Scalar heatFlux_;
    Scalar lambdaSolid_;
    PcKrSwCurve pcKrSwCurve_;
};

} // end namespace Dumux

int main(int argc, char* argv[])
{
    using namespace Dumux;
    using Scalar = HeatPipeReferenceODE::Scalar;
    using Method = Experimental::MultiStage::RungeKuttaExplicitFourthOrder<Scalar>;
    using Variables = HeatPipeReferenceODE::Variables;
    using ODE = HeatPipeReferenceODE;

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv);

    HeatPipeReferenceODE::FluidSystem::init();

    auto ode = std::make_shared<HeatPipeReferenceODE>();
    auto method = std::make_shared<Method>();
    Experimental::ODESolver<HeatPipeReferenceODE> solver(ode, method);

    const auto domainLength = getParam<Scalar>("Reference.DomainLength", 2.4);
    const auto maxStepSize = getParam<Scalar>("Reference.StepSize", 2.5e-4);
    const auto minStepSize = getParam<Scalar>("Reference.MinStepSize", 1e-5);
    std::ofstream output(getParam<std::string>("Reference.OutputFile", "heatpipe_reference.csv"));
    output << std::setprecision(12) << "x,S_liq,p_gas,x_air_gas,T,two_phase\n";

    const auto write = [&](const Scalar x, const Scalar sw, const auto& z, const bool twoPhase)
    {
        output << x << "," << sw << "," << z[ODE::gasPressureIdx] << ","
               << z[ODE::gasAirMoleFractionIdx] << "," << z[ODE::temperatureIdx] << ","
               << twoPhase << "\n";
    };

    // Integrate until the wetting phase dries out (Se -> 0). Steps that would cross
    // the dry-out front or that yield a non-physical state (the explicit scheme can
    // become unstable where the air mole fraction decays to zero) are rejected and
    // repeated with half the step size. This locates the front up to the minimum step size.
    Variables vars(ode->initialState());
    Scalar x = 0.0;
    Scalar stepSize = maxStepSize;
    write(x, ode->saturation(vars.dofs()[ODE::effectiveSaturationIdx]), vars.dofs(), true);
    while (x < domainLength && stepSize >= minStepSize)
    {
        auto trialVars = vars;
        bool stepFailed = false;
        try
        {
            solver.step(trialVars, x, std::min(stepSize, domainLength - x));
        }
        catch (const NumericalProblem&)
        {
            stepFailed = true;
        }

        const auto se = trialVars.dofs()[ODE::effectiveSaturationIdx];
        const auto xa = trialVars.dofs()[ODE::gasAirMoleFractionIdx];
        using std::isfinite;
        if (stepFailed || !isfinite(se) || se <= 1e-6 || !(xa >= 0.0))
        {
            stepSize *= 0.5;
            continue;
        }

        vars = trialVars;
        x = vars.independentVariableLevel().current();
        write(x, ode->saturation(se), vars.dofs(), true);
        stepSize = std::min(2.0*stepSize, maxStepSize);
    }

    // Beyond the dry-out front, the numerical model switches to a gas-only phase
    // state (Sw = 0), and heat is transported by conduction only.
    const auto dryOutPosition = x;
    auto dryState = vars.dofs();
    dryState[ODE::gasAirMoleFractionIdx] = 0.0;
    const auto dryOutTemperature = dryState[ODE::temperatureIdx];
    for (const auto xDry : linspace(dryOutPosition, domainLength, 200))
    {
        dryState[ODE::temperatureIdx] = dryOutTemperature + ode->dryTemperatureGradient()*(xDry - dryOutPosition);
        write(xDry, 0.0, dryState, false);
    }

    std::cout << std::setprecision(6) << "Dry-out front of the semi-analytical solution at x = "
              << dryOutPosition << " m" << std::endl;

    // regression check for the default parameters (params.input)
    using std::abs;
    const auto referenceDryOutPosition = getParam<Scalar>("Problem.ReferenceDryOutPosition");
    if (abs(dryOutPosition - referenceDryOutPosition) > 1e-4)
        DUNE_THROW(Dune::InvalidStateException, "Dry-out front at x = " << dryOutPosition
                    << " m deviates from the expected position x = " << referenceDryOutPosition << " m");

    return 0;
}
