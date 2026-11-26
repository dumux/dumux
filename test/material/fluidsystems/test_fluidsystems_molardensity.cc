// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test that the molar density of a fluid phase is consistent with its
 *        density and composition.
 *
 * The molar density times the average molar mass has to equal the density times
 * the sum of the mole fractions, such that the molar and the mass concentrations
 * of each component agree, also for the mole fractions of a phase that is not
 * present. The models with mole fractions as primary variables balance the
 * components with the molar density and the energy with the density. If the two
 * are inconsistent, a phase transports a different mass in the component balances
 * than in the energy balance.
 */

#include <config.h>

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

#include <dumux/common/parameters.hh>

#include <dumux/material/components/co2.hh>
#include <dumux/material/components/defaultco2table.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/tabulatedcomponent.hh>

#include <dumux/material/fluidstates/compositional.hh>

#include <dumux/material/fluidsystems/brineair.hh>
#include <dumux/material/fluidsystems/brineco2.hh>
#include <dumux/material/fluidsystems/h2oair.hh>
#include <dumux/material/fluidsystems/h2oairmesitylene.hh>
#include <dumux/material/fluidsystems/h2oairxylene.hh>
#include <dumux/material/fluidsystems/h2oheavyoil.hh>
#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/fluidsystems/h2on2o2.hh>

namespace Dumux {

template<class FluidSystem>
using MoleFractions = std::array<double, FluidSystem::numComponents>;

//! Returns 1 if the molar density of the phase is not consistent with its density and composition
template<class FluidSystem>
int checkMolarDensity(const std::string& name, int phaseIdx,
                      double temperature, double pressure,
                      const MoleFractions<FluidSystem>& moleFractions)
{
    CompositionalFluidState<double, FluidSystem> fluidState;
    fluidState.setTemperature(temperature);
    for (int pIdx = 0; pIdx < FluidSystem::numPhases; ++pIdx)
    {
        fluidState.setPressure(pIdx, pressure);
        for (int compIdx = 0; compIdx < FluidSystem::numComponents; ++compIdx)
            fluidState.setMoleFraction(pIdx, compIdx, moleFractions[compIdx]);
    }

    double sumMoleFractions = 0.0;
    for (int compIdx = 0; compIdx < FluidSystem::numComponents; ++compIdx)
        sumMoleFractions += moleFractions[compIdx];

    const double density = FluidSystem::density(fluidState, phaseIdx);
    const double molarDensity = FluidSystem::molarDensity(fluidState, phaseIdx);
    const double averageMolarMass = fluidState.averageMolarMass(phaseIdx);
    const double relativeDifference = std::abs(molarDensity*averageMolarMass - sumMoleFractions*density)
                                      / std::abs(sumMoleFractions*density);

    std::cout << name << ", " << FluidSystem::phaseName(phaseIdx)
              << " at T = " << temperature << " K, p = " << pressure << " Pa, sum of mole fractions "
              << sumMoleFractions << ": relative difference of molarDensity*averageMolarMass and "
              << "density*sumMoleFractions " << relativeDifference;

    // written such that NaN fails
    if (relativeDifference < 1e-12 && molarDensity > 0.0)
    {
        std::cout << std::endl;
        return 0;
    }

    std::cout << " (inconsistent)" << std::endl;
    return 1;
}

//! Liquid and gas phase of the fluid systems of water and one gas, optionally with NaCl in the liquid
template<class FluidSystem>
int checkWaterGasFluidSystem(const std::string& name, bool checkAbsentGas, int gasCompIdx, int naclCompIdx = -1)
{
    FluidSystem::init();

    const int liquidPhaseIdx = FluidSystem::liquidPhaseIdx;
    const int gasPhaseIdx = FluidSystem::gasPhaseIdx;
    const int h2oIdx = FluidSystem::H2OIdx;

    int failures = 0;

    MoleFractions<FluidSystem> liquid{};
    liquid[gasCompIdx] = 1e-5;
    if (naclCompIdx >= 0)
        liquid[naclCompIdx] = 0.05;
    liquid[h2oIdx] = 1.0 - liquid[gasCompIdx] - (naclCompIdx >= 0 ? liquid[naclCompIdx] : 0.0);
    failures += checkMolarDensity<FluidSystem>(name, liquidPhaseIdx, 293.15, 1e5, liquid);

    MoleFractions<FluidSystem> humidAir{};
    humidAir[h2oIdx] = 0.02;
    humidAir[gasCompIdx] = 0.98;
    failures += checkMolarDensity<FluidSystem>(name, gasPhaseIdx, 293.15, 1e5, humidAir);

    // steam-rich gas with the partial pressure of water below its vapor pressure
    MoleFractions<FluidSystem> steam{};
    steam[h2oIdx] = 0.9;
    steam[gasCompIdx] = 0.1;
    failures += checkMolarDensity<FluidSystem>(name, gasPhaseIdx, 523.15, 3.5e6, steam);

    if (checkAbsentGas)
    {
        // gas phase that is not present, with the mole fraction of water given by its vapor pressure
        MoleFractions<FluidSystem> absentGas{};
        absentGas[h2oIdx] = 0.02;
        failures += checkMolarDensity<FluidSystem>(name, gasPhaseIdx, 293.15, 1e5, absentGas);

        // the same during a Newton iteration with a negative gas pressure
        MoleFractions<FluidSystem> negativePressure{};
        negativePressure[h2oIdx] = -0.06;
        failures += checkMolarDensity<FluidSystem>(name, gasPhaseIdx, 290.0, -3e4, negativePressure);
    }

    return failures;
}

//! Water, NAPL and gas phase of the three-phase fluid systems with a pure NAPL phase
template<class FluidSystem>
int checkThreePhaseFluidSystem(const std::string& name, double hotTemperature, double hotPressure)
{
    FluidSystem::init();

    constexpr bool hasAir = FluidSystem::numComponents == 3;
    const int h2oIdx = FluidSystem::H2OIdx;
    const int naplIdx = FluidSystem::NAPLIdx;

    int failures = 0;

    MoleFractions<FluidSystem> water{};
    water[naplIdx] = 1e-5;
    water[h2oIdx] = 1.0 - 1e-5;
    if constexpr (hasAir)
    {
        water[FluidSystem::AirIdx] = 1e-5;
        water[h2oIdx] -= 1e-5;
    }
    failures += checkMolarDensity<FluidSystem>(name, FluidSystem::wPhaseIdx, 293.15, 1e5, water);

    MoleFractions<FluidSystem> napl{};
    napl[naplIdx] = 1.0;
    failures += checkMolarDensity<FluidSystem>(name, FluidSystem::nPhaseIdx, 293.15, 1e5, napl);

    MoleFractions<FluidSystem> humidGas{};
    humidGas[h2oIdx] = 0.02;
    if constexpr (hasAir)
    {
        humidGas[naplIdx] = 1e-3;
        humidGas[FluidSystem::AirIdx] = 1.0 - 0.02 - 1e-3;
    }
    else
        humidGas[naplIdx] = 0.98;
    failures += checkMolarDensity<FluidSystem>(name, FluidSystem::gPhaseIdx, 293.15, 1e5, humidGas);

    MoleFractions<FluidSystem> steam{};
    steam[h2oIdx] = 0.9;
    if constexpr (hasAir)
    {
        steam[naplIdx] = 0.01;
        steam[FluidSystem::AirIdx] = 0.09;
    }
    else
        steam[naplIdx] = 0.1;
    failures += checkMolarDensity<FluidSystem>(name, FluidSystem::gPhaseIdx, hotTemperature, hotPressure, steam);

    // gas phase that is not present, with the mole fractions given by the vapor pressures
    MoleFractions<FluidSystem> absentGas{};
    absentGas[h2oIdx] = 0.02;
    absentGas[naplIdx] = 1e-3;
    failures += checkMolarDensity<FluidSystem>(name, FluidSystem::gPhaseIdx, 293.15, 1e5, absentGas);

    // the same during a Newton iteration with a negative gas pressure
    MoleFractions<FluidSystem> negativePressure{};
    negativePressure[h2oIdx] = -0.06;
    failures += checkMolarDensity<FluidSystem>(name, FluidSystem::gPhaseIdx, 290.0, -3e4, negativePressure);

    return failures;
}

//! Liquid and gas phase of the brine-CO2 fluid systems, within the range of the default CO2 tables
template<class FluidSystem>
int checkBrineCO2FluidSystem(const std::string& name)
{
    FluidSystem::init();

    constexpr bool hasNaCl = FluidSystem::numComponents == 3;
    const int h2oIdx = FluidSystem::comp0Idx;
    const int co2Idx = FluidSystem::CO2Idx;

    int failures = 0;

    MoleFractions<FluidSystem> liquid{};
    liquid[co2Idx] = 0.02;
    liquid[h2oIdx] = 0.98;
    if constexpr (hasNaCl)
    {
        liquid[FluidSystem::comp2Idx] = 0.05;
        liquid[h2oIdx] -= 0.05;
    }
    failures += checkMolarDensity<FluidSystem>(name, FluidSystem::liquidPhaseIdx, 313.15, 1e7, liquid);

    // supercritical CO2 with the partial pressure of water below its vapor pressure
    MoleFractions<FluidSystem> supercritical{};
    supercritical[h2oIdx] = 5e-4;
    supercritical[co2Idx] = 1.0 - 5e-4;
    failures += checkMolarDensity<FluidSystem>(name, FluidSystem::gasPhaseIdx, 313.15, 1e7, supercritical);

    // gas with the partial pressure of water below its vapor pressure
    MoleFractions<FluidSystem> humidGas{};
    humidGas[h2oIdx] = 0.02;
    humidGas[co2Idx] = 0.98;
    failures += checkMolarDensity<FluidSystem>(name, FluidSystem::gasPhaseIdx, 338.15, 1e6, humidGas);

    // gas phase that is not present
    MoleFractions<FluidSystem> absentGas{};
    absentGas[h2oIdx] = 5e-4;
    absentGas[co2Idx] = 0.5;
    failures += checkMolarDensity<FluidSystem>(name, FluidSystem::gasPhaseIdx, 313.15, 1e7, absentGas);

    return failures;
}

} // end namespace Dumux

int main()
{
    using namespace Dumux;

    using Scalar = double;
    using H2O = Components::H2O<Scalar>;
    using SimpleH2O = Components::SimpleH2O<Scalar>;
    using TabulatedH2O = Components::TabulatedComponent<H2O>;
    using CO2 = Components::CO2<Scalar, GeneratedCO2Tables::CO2Tables>;

    Parameters::init([](auto& params){ params["Brine.Salinity"] = "0.1"; });
    std::cout << std::scientific << std::setprecision(3);

    int failures = 0;

    // H2O -- Air
    {   using FS = FluidSystems::H2OAir<Scalar, H2O>;
        failures += checkWaterGasFluidSystem<FS>("H2OAir<H2O>", true, FS::AirIdx); }
    {   using FS = FluidSystems::H2OAir<Scalar, H2O, FluidSystems::H2OAirDefaultPolicy</*fastButSimplifiedRelations=*/true>>;
        failures += checkWaterGasFluidSystem<FS>("H2OAir<H2O, simplified>", false, FS::AirIdx); }
    {   using FS = FluidSystems::H2OAir<Scalar, TabulatedH2O>;
        failures += checkWaterGasFluidSystem<FS>("H2OAir<TabulatedH2O>", true, FS::AirIdx); }
    {   using FS = FluidSystems::H2OAir<Scalar, SimpleH2O>;
        failures += checkWaterGasFluidSystem<FS>("H2OAir<SimpleH2O>", true, FS::AirIdx); }

    // Brine -- Air
    {   using FS = FluidSystems::BrineAir<Scalar, H2O>;
        failures += checkWaterGasFluidSystem<FS>("BrineAir<H2O>", true, FS::AirIdx, FS::NaClIdx); }
    {   using FS = FluidSystems::BrineAir<Scalar, H2O, FluidSystems::BrineAirDefaultPolicy</*fastButSimplifiedRelations=*/true>>;
        failures += checkWaterGasFluidSystem<FS>("BrineAir<H2O, simplified>", false, FS::AirIdx, FS::NaClIdx); }
    {   using FS = FluidSystems::BrineAir<Scalar, SimpleH2O>;
        failures += checkWaterGasFluidSystem<FS>("BrineAir<SimpleH2O>", true, FS::AirIdx, FS::NaClIdx); }

    // H2O -- N2
    {   using FS = FluidSystems::H2ON2<Scalar>;
        failures += checkWaterGasFluidSystem<FS>("H2ON2", true, FS::N2Idx); }
    {   using FS = FluidSystems::H2ON2<Scalar, FluidSystems::H2ON2DefaultPolicy</*fastButSimplifiedRelations=*/true>>;
        failures += checkWaterGasFluidSystem<FS>("H2ON2<simplified>", false, FS::N2Idx); }

    // H2O -- N2 -- O2
    {   using FS = FluidSystems::H2ON2O2<Scalar>;
        failures += checkWaterGasFluidSystem<FS>("H2ON2O2", true, FS::N2Idx); }
    {   using FS = FluidSystems::H2ON2O2<Scalar, FluidSystems::H2ON2O2DefaultPolicy</*fastButSimplifiedRelations=*/true>>;
        failures += checkWaterGasFluidSystem<FS>("H2ON2O2<simplified>", false, FS::N2Idx); }

    // three phases
    {   using FS = FluidSystems::H2OAirMesitylene<Scalar>;
        failures += checkThreePhaseFluidSystem<FS>("H2OAirMesitylene", 373.15, 1e5); }
    {   using FS = FluidSystems::H2OAirXylene<Scalar>;
        failures += checkThreePhaseFluidSystem<FS>("H2OAirXylene", 373.15, 1e5); }
    {   using FS = FluidSystems::H2OHeavyOil<Scalar>;
        failures += checkThreePhaseFluidSystem<FS>("H2OHeavyOil", 523.15, 3.5e6); }

    // Brine -- CO2
    {   using FS = FluidSystems::BrineCO2<Scalar, CO2, TabulatedH2O, FluidSystems::BrineCO2DefaultPolicy</*constantSalinity=*/true>>;
        failures += checkBrineCO2FluidSystem<FS>("BrineCO2<TabulatedH2O, constant salinity>"); }
    {   using FS = FluidSystems::BrineCO2<Scalar, CO2, TabulatedH2O, FluidSystems::BrineCO2DefaultPolicy</*constantSalinity=*/false>>;
        failures += checkBrineCO2FluidSystem<FS>("BrineCO2<TabulatedH2O>"); }
    {   using FS = FluidSystems::BrineCO2<Scalar, CO2, TabulatedH2O, FluidSystems::BrineCO2DefaultPolicy</*constantSalinity=*/false, /*fastButSimplifiedRelations=*/true>>;
        failures += checkBrineCO2FluidSystem<FS>("BrineCO2<TabulatedH2O, simplified>"); }
    {   using FS = FluidSystems::BrineCO2<Scalar, CO2, H2O, FluidSystems::BrineCO2DefaultPolicy</*constantSalinity=*/false>>;
        failures += checkBrineCO2FluidSystem<FS>("BrineCO2<H2O>"); }
    {   using FS = FluidSystems::BrineCO2<Scalar, CO2, SimpleH2O, FluidSystems::BrineCO2DefaultPolicy</*constantSalinity=*/false>>;
        failures += checkBrineCO2FluidSystem<FS>("BrineCO2<SimpleH2O>"); }

    if (failures > 0)
    {
        std::cout << failures << " inconsistent molar densities" << std::endl;
        return 1;
    }

    return 0;
}
