// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup RichardsTests
 * \brief Checks that the extended Richards model's equilibrium vapor
 *        pressure follows Kelvin's equation when the fluid system provides
 *        a fluid-state dependent vapor pressure, and falls back to the
 *        plain component vapor pressure otherwise.
 */
#include <config.h>

#include <cmath>
#include <iostream>
#include <string>

#include <dumux/material/constants.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/fluidsystems/h2oair.hh>
#include <dumux/material/fluidstates/immiscible.hh>
#include <dumux/porousmediumflow/richardsextended/primaryvariableswitch.hh>

namespace Dumux::Test {

// a minimal fluid system that does not offer a fluid-state dependent vapor pressure
template<class Scalar>
struct WaterOnlyFluidSystem
{
    using H2O = Dumux::Components::H2O<Scalar>;
    static constexpr int liquidPhaseIdx = 0;
};

template<class Scalar>
bool fuzzyEqual(Scalar a, Scalar b, Scalar tolerance = 1e-12)
{
    using std::abs;
    return abs(a-b) <= tolerance*abs(b);
}

} // end namespace Dumux::Test

int main()
{
    using namespace Dumux;

    using Scalar = double;
    using H2O = Components::H2O<Scalar>;
    using H2OAirNoKelvin = FluidSystems::H2OAir<Scalar, H2O, FluidSystems::H2OAirDefaultPolicy<>, /*useKelvinVaporPressure=*/false>;
    using H2OAirKelvin = FluidSystems::H2OAir<Scalar, H2O, FluidSystems::H2OAirDefaultPolicy<>, /*useKelvinVaporPressure=*/true>;

    const Scalar temperature = 293.15;
    const Scalar liquidPressure = 1e5;
    const Scalar pc = 5e4;
    const Scalar gasPressure = liquidPressure + pc;

    int failures = 0;

    // the fallback path is used for a fluid system that offers no fluid-state
    // dependent vapor pressure: the equilibrium vapor pressure is that of pure water
    {
        using FluidSystem = Test::WaterOnlyFluidSystem<Scalar>;
        using FluidState = ImmiscibleFluidState<Scalar, H2OAirNoKelvin>;

        FluidState fs;
        fs.setTemperature(temperature);
        fs.setPressure(H2OAirNoKelvin::liquidPhaseIdx, liquidPressure);
        fs.setPressure(H2OAirNoKelvin::gasPhaseIdx, gasPressure);
        fs.setWettingPhase(H2OAirNoKelvin::liquidPhaseIdx);

        const Scalar actual = Detail::ExtendedRichards::equilibriumVaporPressure<FluidSystem>(fs);
        const Scalar expected = H2O::vaporPressure(temperature);

        if (!Test::fuzzyEqual(actual, expected))
        {
            std::cerr << "Fallback vapor pressure wrong: expected " << expected << ", got " << actual << std::endl;
            ++failures;
        }
    }

    // H2OAir without Kelvin's equation provides a fluid-state dependent vapor pressure
    // that must still equal the pure water vapor pressure (no capillary pressure effect)
    {
        using FluidSystem = H2OAirNoKelvin;
        using FluidState = ImmiscibleFluidState<Scalar, FluidSystem>;

        FluidState fs;
        fs.setTemperature(temperature);
        fs.setPressure(FluidSystem::liquidPhaseIdx, liquidPressure);
        fs.setPressure(FluidSystem::gasPhaseIdx, gasPressure);
        fs.setWettingPhase(FluidSystem::liquidPhaseIdx);

        const Scalar actual = Detail::ExtendedRichards::equilibriumVaporPressure<FluidSystem>(fs);
        const Scalar expected = H2O::vaporPressure(temperature);

        if (!Test::fuzzyEqual(actual, expected))
        {
            std::cerr << "H2OAir (no Kelvin) vapor pressure wrong: expected " << expected << ", got " << actual << std::endl;
            ++failures;
        }
    }

    // H2OAir with Kelvin's equation lowers the vapor pressure according to
    // p_v = p_sat(T) * exp(-pc*M/(rho_w*R*T))
    {
        using FluidSystem = H2OAirKelvin;
        using FluidState = ImmiscibleFluidState<Scalar, FluidSystem>;

        FluidState fs;
        fs.setTemperature(temperature);
        fs.setPressure(FluidSystem::liquidPhaseIdx, liquidPressure);
        fs.setPressure(FluidSystem::gasPhaseIdx, gasPressure);
        fs.setWettingPhase(FluidSystem::liquidPhaseIdx);

        const Scalar actual = Detail::ExtendedRichards::equilibriumVaporPressure<FluidSystem>(fs);

        using std::exp;
        const Scalar molarMass = H2O::molarMass();
        const Scalar liquidDensity = H2O::liquidDensity(temperature, liquidPressure);
        const Scalar expected = H2O::vaporPressure(temperature)
                                 * exp(-pc*molarMass/(liquidDensity*Constants<Scalar>::R*temperature));

        if (!Test::fuzzyEqual(actual, expected))
        {
            std::cerr << "H2OAir (Kelvin) vapor pressure wrong: expected " << expected << ", got " << actual << std::endl;
            ++failures;
        }

        // Kelvin's equation must lower the vapor pressure below the flat-interface value
        if (!(actual < H2O::vaporPressure(temperature)))
        {
            std::cerr << "Kelvin vapor pressure is not below the flat-interface vapor pressure" << std::endl;
            ++failures;
        }
    }

    return failures;
}
