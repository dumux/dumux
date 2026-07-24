// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/experimental/ode/odesolver.hh>

namespace Dumux {

class HeatPipeReferenceODE
{
public:
    using Scalar = double;
    using SolutionVector = Dune::FieldVector<Scalar, 4>;
    using ResidualType = SolutionVector;
    using JacobianMatrix = Dune::FieldMatrix<Scalar, 4, 4>;
    using Variables = Experimental::Variables<SolutionVector>;

    enum Component
    {
        effectiveSaturationIdx,
        gasPressureIdx,
        gasAirMoleFractionIdx,
        temperatureIdx
    };

    SolutionVector initialState() const
    {
        SolutionVector result;
        result[effectiveSaturationIdx] = (swBc_ - swr_)/(1.0 - swr_);
        result[gasPressureIdx] = pgBc_;
        result[gasAirMoleFractionIdx] = gasAirMoleFractionBoundary_(result[effectiveSaturationIdx]);
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

        const auto sw = swr_ + se*(1.0 - swr_);
        const auto sg = 1.0 - sw;

        const auto dPm = millingtonQuirk_(sg, binaryDiffusionCoefficient_(temperature, pg));
        const auto lambda = lambdaDry_ + std::sqrt(sw)*(lambdaWet_ - lambdaDry_);
        const auto pc = capillaryPressure_(se);
        const auto dpcDse = capillaryPressureDerivative_(se);

        const auto krg = gasRelativePermeability_(se);
        const auto krl = liquidRelativePermeability_(se);
        const auto rhoGa = molarMassAir_*xa*pg/(gasConstant_*temperature);
        const auto rhoGw = molarMassWater_*(1.0 - xa)*pg/(gasConstant_*temperature);
        const auto rhoG = rhoGa + rhoGw;
        const auto muG = xa*muAir_ + (1.0 - xa)*muSteam_;
        const auto nuG = muG/rhoG;
        const auto nuW = muWater_/rhoWater_;
        const auto beta = nuW/nuG;
        const auto alpha = 1.0 + pc/(rhoWater_*latentHeat_);
        const auto xi = (1.0/krg)*(1.0 + rhoWater_*gasConstant_*temperature
                                                /(pg*molarMassWater_*(1.0 - xa)))
                        + beta/krl;
        const auto delta = rhoWater_*latentHeat_*latentHeat_*permeability_*alpha/(lambda*nuG*temperature);
        const auto zeta = permeability_*rhoWater_*gasConstant_*temperature/(molarMassWater_*rhoG*nuG*dPm)
                          *xa/(1.0 - xa)
                          *(pg*molarMassWater_/(rhoWater_*gasConstant_*temperature) + 1.0/(1.0 - xa));
        const auto eta = delta/(delta + xi + zeta);

        rhs[effectiveSaturationIdx] = -(1.0/(1.0 - xa) + beta*krg/krl)
                                      *eta*heatFlux_*nuG/(permeability_*latentHeat_*krg*dpcDse);
        rhs[gasPressureIdx] = -(eta*heatFlux_*nuG/(permeability_*latentHeat_*krg))/(1.0 - xa);
        rhs[gasAirMoleFractionIdx] = eta*heatFlux_*xa/(latentHeat_*dPm*rhoG*(1.0 - xa));
        rhs[temperatureIdx] = -heatFlux_*(1.0 - eta)/lambda;
    }

private:
    Scalar p0Leverett_() const
    { return std::sqrt(porosity_/permeability_); }

    Scalar capillaryPressure_(const Scalar se) const
    {
        const auto oneMinusSe = 1.0 - se;
        return p0Leverett_()*surfaceTension_*(1.417*oneMinusSe
                                             - 2.120*oneMinusSe*oneMinusSe
                                             + 1.263*oneMinusSe*oneMinusSe*oneMinusSe);
    }

    Scalar capillaryPressureDerivative_(const Scalar se) const
    {
        const auto eps = 1e-8;
        return (capillaryPressure_(se + eps) - capillaryPressure_(se - eps))/(2.0*eps);
    }

    Scalar liquidRelativePermeability_(const Scalar se) const
    { return se*se*se; }

    Scalar gasRelativePermeability_(const Scalar se) const
    {
        const auto oneMinusSe = 1.0 - se;
        return oneMinusSe*oneMinusSe*oneMinusSe;
    }

    Scalar binaryDiffusionCoefficient_(const Scalar temperature, const Scalar pressure) const
    {
        return diffusionReference_*diffusionPressureReference_/pressure
               *std::pow(temperature/diffusionTemperatureReference_, diffusionExponent_);
    }

    Scalar millingtonQuirk_(const Scalar sg, const Scalar binaryDiffusionCoefficient) const
    {
        const auto positiveSg = std::max(sg, 0.0);
        return porosity_*std::pow(positiveSg, 3)
               *std::pow(porosity_*positiveSg, 1.0/3.0)*binaryDiffusionCoefficient;
    }

    Scalar gasAirMoleFractionBoundary_(const Scalar se) const
    {
        const auto pc = capillaryPressure_(se);
        return 1.0 - pressureReference_/pgBc_
                     *std::exp((1.0/temperatureReference_ - 1.0/tBc_)*latentHeat_*molarMassWater_/gasConstant_
                               - pc*molarMassWater_/(rhoWater_*gasConstant_*tBc_));
    }

    static constexpr Scalar permeability_ = 1.0e-12;
    static constexpr Scalar porosity_ = 0.4;
    static constexpr Scalar swr_ = 0.15;
    static constexpr Scalar rhoWater_ = 958.4;
    static constexpr Scalar muWater_ = 2.938e-4;
    static constexpr Scalar muAir_ = 2.08e-5;
    static constexpr Scalar muSteam_ = 1.2e-5;
    static constexpr Scalar surfaceTension_ = 0.05878;
    static constexpr Scalar molarMassWater_ = 0.018016;
    static constexpr Scalar molarMassAir_ = 0.02897;
    static constexpr Scalar gasConstant_ = 8.3144621;
    static constexpr Scalar latentHeat_ = 2.258e6;
    static constexpr Scalar pressureReference_ = 101325.0;
    static constexpr Scalar temperatureReference_ = 373.15;
    static constexpr Scalar lambdaWet_ = 1.586;
    static constexpr Scalar lambdaDry_ = 0.4278;
    static constexpr Scalar diffusionExponent_ = 1.8;
    static constexpr Scalar diffusionReference_ = 2.13e-5;
    static constexpr Scalar diffusionPressureReference_ = 1.0e5;
    static constexpr Scalar diffusionTemperatureReference_ = 273.15;
    static constexpr Scalar pgBc_ = 1.013e5;
    static constexpr Scalar swBc_ = 0.99;
    static constexpr Scalar tBc_ = 341.75;
    static constexpr Scalar heatFlux_ = -100.0;
};

} // end namespace Dumux

int main(int argc, char* argv[])
{
    using namespace Dumux;
    using Scalar = double;
    using Method = Experimental::MultiStage::RungeKuttaExplicitFourthOrder<Scalar>;
    using SolutionVector = HeatPipeReferenceODE::SolutionVector;

    Dumux::initialize(argc, argv);
    Dumux::Parameters::init(argc, argv);

    const auto expectNear = [] (const auto value,
                                const auto reference,
                                const auto tolerance,
                                const std::string& message)
    {
        using std::abs;
        if (abs(value - reference) > tolerance)
            DUNE_THROW(Dune::InvalidStateException,
                       message << " (value: " << value << ", reference: " << reference << ")");
    };

    auto ode = std::make_shared<HeatPipeReferenceODE>();
    auto method = std::make_shared<Method>();
    Experimental::ODESolver<HeatPipeReferenceODE> solver(ode, method);

    Experimental::Variables<SolutionVector> vars(ode->initialState());
    solver.solve(vars, 0.25, 2.5e-4);

    expectNear(vars.dofs()[HeatPipeReferenceODE::effectiveSaturationIdx], 0.98803128090343983, 1e-8,
               "Heatpipe ODE solve failed for the effective saturation at x=0.25 m");
    expectNear(vars.dofs()[HeatPipeReferenceODE::gasPressureIdx], 101310.36945347066, 1e-5,
               "Heatpipe ODE solve failed for the gas pressure at x=0.25 m");
    expectNear(vars.dofs()[HeatPipeReferenceODE::gasAirMoleFractionIdx], 0.43520555395030663, 1e-8,
               "Heatpipe ODE solve failed for the gas-phase air mole fraction at x=0.25 m");
    expectNear(vars.dofs()[HeatPipeReferenceODE::temperatureIdx], 357.57126217848389, 1e-8,
               "Heatpipe ODE solve failed for the temperature at x=0.25 m");

    solver.solve(vars, 1.0, 2.5e-4);

    expectNear(vars.dofs()[HeatPipeReferenceODE::effectiveSaturationIdx], 0.38566288370571478, 1e-6,
               "Heatpipe ODE solve failed for the effective saturation at x=1.0 m");
    expectNear(vars.dofs()[HeatPipeReferenceODE::gasPressureIdx], 114144.93970857002, 1e-1,
               "Heatpipe ODE solve failed for the gas pressure at x=1.0 m");
    expectNear(vars.dofs()[HeatPipeReferenceODE::temperatureIdx], 376.58683441201634, 5e-5,
               "Heatpipe ODE solve failed for the temperature at x=1.0 m");
    expectNear(vars.dofs()[HeatPipeReferenceODE::gasAirMoleFractionIdx], 0.0, 1e-8,
               "Heatpipe ODE solve failed for the depleted gas-phase air mole fraction at x=1.0 m");

    std::cout << "Heatpipe ODE solver test passed" << std::endl;

    return 0;
}
