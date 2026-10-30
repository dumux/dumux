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
 * towards the heat source until the wetting phase dries out. The material
 * properties, the extensions to an ideal gas mixture with a mole-fraction
 * weighted viscosity and the boundary state follow the supplementary code of
 * Huang et al. (2015). They are implemented here independently of the numerical
 * model, which reproduces them with referencefluidsystem.hh,
 * thermalconductivity.hh and diffusivity.hh.
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
    {}

    SolutionVector initialState() const
    {
        SolutionVector result;
        result[effectiveSaturationIdx] = (swBc_ - swr_)/(1.0 - swr_);
        result[gasPressureIdx] = pgBc_;
        result[gasAirMoleFractionIdx] = xaBc_;
        result[temperatureIdx] = boundaryTemperature();
        return result;
    }

    //! The temperature in local equilibrium with the boundary state, eq. (11) in Huang et al. (2015)
    Scalar boundaryTemperature() const
    {
        const auto se = (swBc_ - swr_)/(1.0 - swr_);
        const auto pc = capillaryPressure_(se);
        return boilingTemperature_*(1.0 + (pc - xaBc_*pgBc_)/(latentHeat_*liquidDensity_))
               /(1.0 - boilingTemperature_*gasConstant_/(latentHeat_*molarMassWater_)
                       *std::log(pgBc_*(1.0 - xaBc_)/boilingPressure_));
    }

    void rhs(const Variables& vars, ResidualType& rhs) const
    {
        const auto& z = vars.dofs();
        const auto se = z[effectiveSaturationIdx];
        const auto pg = z[gasPressureIdx];
        const auto xa = z[gasAirMoleFractionIdx];
        const auto temperature = z[temperatureIdx];

        const auto sw = saturation(se);
        const auto pc = capillaryPressure_(se);
        const auto dpcDse = capillaryPressureDerivative_(se);
        const auto krg = (1.0 - se)*(1.0 - se)*(1.0 - se);
        const auto krl = se*se*se;
        const auto dPm = diffusionCoefficient_;
        const auto lambda = thermalConductivity_(sw);

        const auto rhoG = (molarMassAir_*xa + molarMassWater_*(1.0 - xa))*pg/(gasConstant_*temperature);
        const auto muG = xa*airViscosity_ + (1.0 - xa)*vaporViscosity_;
        const auto nuG = muG/rhoG;
        const auto nuW = liquidViscosity_/liquidDensity_;
        const auto beta = nuW/nuG;
        const auto alpha = 1.0 + (pc - xa*pg)/(liquidDensity_*latentHeat_);
        const auto xi = (1.0/krg)*(1.0 + liquidDensity_*gasConstant_*temperature
                                         /(pg*molarMassWater_*(1.0 - xa)))
                        + beta/krl;
        const auto delta = liquidDensity_*latentHeat_*latentHeat_*permeability_*alpha/(lambda*nuG*temperature);
        const auto zeta = permeability_*liquidDensity_*gasConstant_*temperature/(molarMassWater_*rhoG*nuG*dPm)
                          *xa/(1.0 - xa)
                          *(pg*molarMassWater_/(liquidDensity_*gasConstant_*temperature) + 1.0/(1.0 - xa));
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
    { return -heatFlux_/thermalConductivity_(0.0); }

private:
    //! Leverett function as a function of the effective saturation
    Scalar capillaryPressure_(const Scalar se) const
    {
        const auto sn = 1.0 - se;
        return surfaceTension_*std::sqrt(porosity_/permeability_)*((1.263*sn - 2.120)*sn + 1.417)*sn;
    }

    Scalar capillaryPressureDerivative_(const Scalar se) const
    {
        const auto sn = 1.0 - se;
        return -surfaceTension_*std::sqrt(porosity_/permeability_)*((3.0*1.263*sn - 2.0*2.120)*sn + 1.417);
    }

    Scalar thermalConductivity_(const Scalar sw) const
    { return lambdaDry_ + std::sqrt(std::max(sw, 0.0))*(lambdaWet_ - lambdaDry_); }

    // material properties of Huang et al. (2015), supplementary code
    static constexpr Scalar porosity_ = 0.35;
    static constexpr Scalar swr_ = 0.0;
    static constexpr Scalar surfaceTension_ = 0.05878; // N/m
    static constexpr Scalar lambdaWet_ = 1.13; // W/(m K)
    static constexpr Scalar lambdaDry_ = 0.582; // W/(m K)
    static constexpr Scalar liquidDensity_ = 1000.0; // kg/m^3
    static constexpr Scalar liquidViscosity_ = 2.938e-4; // Pa s
    static constexpr Scalar airViscosity_ = 2.094e-5; // Pa s
    static constexpr Scalar vaporViscosity_ = 1.2e-5; // Pa s
    static constexpr Scalar diffusionCoefficient_ = 2.6e-6; // m^2/s, macroscopic (porous medium) value
    static constexpr Scalar latentHeat_ = 2.258e6; // J/kg
    static constexpr Scalar boilingTemperature_ = 373.15; // K
    static constexpr Scalar boilingPressure_ = 101325.0; // Pa
    static constexpr Scalar gasConstant_ = 8.314; // J/(mol K)
    static constexpr Scalar molarMassWater_ = 0.018016; // kg/mol
    static constexpr Scalar molarMassAir_ = 0.028; // kg/mol
    // state at the left boundary, see problem.hh
    static constexpr Scalar pgBc_ = 101325.0;
    static constexpr Scalar swBc_ = 0.99;
    static constexpr Scalar xaBc_ = 0.7;

    Scalar permeability_;
    Scalar heatFlux_;
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

    auto ode = std::make_shared<HeatPipeReferenceODE>();

    // the numerical model prescribes the same boundary temperature (params.input)
    using std::abs;
    const auto leftTemperature = getParam<Scalar>("Problem.LeftTemperature");
    if (abs(ode->boundaryTemperature() - leftTemperature) > 1e-3)
        DUNE_THROW(Dune::InvalidStateException, "Boundary temperature " << ode->boundaryTemperature()
                    << " K of the reference deviates from Problem.LeftTemperature = " << leftTemperature << " K");
    auto method = std::make_shared<Method>();
    Experimental::ODESolver<HeatPipeReferenceODE> solver(ode, method);

    const auto domainLength = getParam<Scalar>("Reference.DomainLength", 2.0);
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
    // become unstable where the air mole fraction decays to zero, and it can overshoot
    // to a larger saturation close to the front, where dSe/dx diverges) are rejected and
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
        if (stepFailed || !isfinite(se) || se <= 1e-6 || !(xa >= 0.0)
            || se >= vars.dofs()[ODE::effectiveSaturationIdx])
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
    const auto referenceDryOutPosition = getParam<Scalar>("Problem.ReferenceDryOutPosition");
    if (abs(dryOutPosition - referenceDryOutPosition) > 1e-4)
        DUNE_THROW(Dune::InvalidStateException, "Dry-out front at x = " << dryOutPosition
                    << " m deviates from the expected position x = " << referenceDryOutPosition << " m");

    return 0;
}
