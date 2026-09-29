// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterTests
 * \brief Green-Ampt infiltration against its closed-form solution.
 */
#include <config.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

#include <dumux/freeflow/shallowwater/surfacerunoff/greenampt.hh>

namespace {

int failures = 0;

void check(const std::string& what, const double value, const double reference, const double tolerance)
{
    using std::abs, std::max;
    const auto error = abs(value - reference)/max(1e-30, abs(reference));
    if (!(error <= tolerance))
    {
        std::cerr << "FAIL " << what << ": got " << std::setprecision(10) << value
                  << ", expected " << reference << " (relative error " << error
                  << " > " << tolerance << ")\n";
        ++failures;
    }
    else
        std::cout << "  ok  " << what << " = " << std::setprecision(8) << value << "\n";
}

void checkTrue(const std::string& what, const bool condition)
{
    if (!condition)
    {
        std::cerr << "FAIL " << what << "\n";
        ++failures;
    }
    else
        std::cout << "  ok  " << what << "\n";
}

/*!
 * \brief Independent solve of the ponded Green-Ampt equation by bisection.
 *
 * Deliberately not the Newton iteration the model uses, so agreement means the two
 * derivations agree and not merely that one implementation is self-consistent.
 */
double cumulativeInfiltration(const double conductivity, const double suctionTimesDeficit,
                              const double time)
{
    const auto residual = [&](const double f)
    { return f - suctionTimesDeficit*std::log1p(f/suctionTimesDeficit) - conductivity*time; };

    double low = 0.0, high = conductivity*time + 2.0*suctionTimesDeficit + 1.0;
    while (residual(high) < 0.0)
        high *= 2.0;
    for (int i = 0; i < 200; ++i)
    {
        const auto mid = 0.5*(low + high);
        (residual(mid) < 0.0 ? low : high) = mid;
    }
    return 0.5*(low + high);
}

} // end anonymous namespace

int main()
{
    using Soil = Dumux::SurfaceRunoff::GreenAmptSoil<double>;

    // a silt loam: Ks 20 mm/h, suction 0.17 m, deficit 0.25, 1 m of soil
    const double conductivity = 20.0e-3/3600.0;
    const double suction = 0.17, deficit = 0.25, soilDepth = 1.0;
    const double suctionTimesDeficit = suction*deficit;
    const Soil soil(conductivity, suction, deficit, soilDepth);

    std::cout << "-- store\n";
    check("capacity", soil.capacity(), soilDepth*deficit, 1e-14);

    std::cout << "\n-- rainfall at or below Ks never ponds\n";
    const auto gentle = 0.5*conductivity;
    checkTrue("time to ponding is infinite",
              std::isinf(soil.timeToPonding(gentle)));
    check("all of it infiltrates", soil.infiltration(0.0, gentle, 3600.0), gentle*3600.0, 1e-14);

    std::cout << "\n-- ponding time is Fp/i, the pre-ponding phase being rainfall-limited\n";
    const auto intensity = 60.0e-3/3600.0; // 60 mm/h, three times Ks
    const auto pondingDepth = suctionTimesDeficit*conductivity/(intensity - conductivity);
    check("wetted depth at ponding", soil.wettedDepthAtPonding(intensity), pondingDepth, 1e-12);
    check("time to ponding", soil.timeToPonding(intensity), pondingDepth/intensity, 1e-12);
    check("capacity has fallen to the rainfall rate there",
          soil.capacityRate(pondingDepth), intensity, 1e-12);
    // up to that moment nothing runs off
    const auto justBefore = 0.99*soil.timeToPonding(intensity);
    check("rainfall-limited before ponding",
          soil.infiltration(0.0, intensity, justBefore), intensity*justBefore, 1e-12);

    std::cout << "\n-- ponded infiltration matches the closed form\n";
    // supply far above Ks so that ponding is effectively immediate
    const Soil deep(conductivity, suction, deficit, 1e6);
    const auto flooded = 1e4*conductivity;
    for (const double hours : {0.5, 2.0, 12.0})
    {
        const auto dt = hours*3600.0;
        const auto pondedFrom = deep.wettedDepthAtPonding(flooded);
        const auto lead = pondedFrom/flooded;
        const auto modelled = deep.infiltration(0.0, flooded, dt);
        // the closed form starts at the instant of ponding, the model a moment earlier
        const auto reference = cumulativeInfiltration(conductivity, suctionTimesDeficit, dt - lead);
        check("F after " + std::to_string(int(hours)) + " h", modelled, reference, 1e-6);
    }

    std::cout << "\n-- marching in small steps agrees with one big step\n";
    for (const int steps : {1, 10, 1000})
    {
        const auto total = 6.0*3600.0;
        const auto dt = total/steps;
        double wetted = 0.0;
        for (int i = 0; i < steps; ++i)
            wetted += deep.infiltration(wetted, intensity, dt);
        // exact for constant supply: the algorithm splits the step at the ponding instant
        const auto reference = deep.infiltration(0.0, intensity, total);
        check(std::to_string(steps) + " steps", wetted, reference, 1e-9);
    }

    std::cout << "\n-- the store fills and then everything runs off\n";
    {
        double wetted = 0.0, infiltrated = 0.0, rain = 0.0;
        bool withinSupply = true;
        const auto dt = 600.0;
        for (int i = 0; i < 500; ++i)
        {
            const auto entered = soil.infiltration(wetted, intensity, dt);
            withinSupply = withinSupply && entered <= intensity*dt + 1e-15;
            wetted += entered;
            infiltrated += entered;
            rain += intensity*dt;
        }
        checkTrue("never takes more than the supply", withinSupply);
        check("wetted depth stops at the capacity", wetted, soil.capacity(), 1e-12);
        checkTrue("the full store took all it could", infiltrated <= soil.capacity() + 1e-12);
        check("runoff is the rest", rain - infiltrated, rain - soil.capacity(), 1e-12);
        check("no further infiltration once full", soil.infiltration(wetted, intensity, dt), 0.0, 1e-30);
    }

    std::cout << "\n-- daily rainfall intensity does not pond a real soil\n";
    {
        // the wettest Hans day, 114.8 mm spread over 24 h
        const auto daily = 114.8e-3/86400.0;
        checkTrue("114.8 mm/day is below Ks", daily < conductivity);
        checkTrue("so it never ponds", std::isinf(soil.timeToPonding(daily)));
        // and yet the store still fills, which is what generates the runoff
        double wetted = 0.0;
        for (int day = 0; day < 11; ++day)
            wetted += soil.infiltration(wetted, daily, 86400.0);
        checkTrue("11 such days still fill the store", wetted >= soil.capacity() - 1e-12);
    }

    std::cout << "\n-- degenerate parameters\n";
    {
        const Soil noSuction(conductivity, 0.0, deficit, soilDepth);
        check("without suction the capacity is Ks",
              noSuction.infiltration(0.0, intensity, 3600.0), conductivity*3600.0, 1e-14);
        const Soil noStore(conductivity, suction, deficit, 0.0);
        check("a zero store takes nothing", noStore.infiltration(0.0, intensity, 3600.0), 0.0, 1e-30);
        check("zero supply takes nothing", soil.infiltration(0.0, 0.0, 3600.0), 0.0, 1e-30);
        const Soil sealed(0.0, suction, deficit, soilDepth);
        check("a sealed surface takes nothing", sealed.infiltration(0.0, intensity, 3600.0), 0.0, 1e-30);
        check("a sealed surface takes nothing once wetted", sealed.infiltration(0.01, intensity, 3600.0), 0.0, 1e-30);
    }

    if (failures > 0)
    {
        std::cerr << "\n" << failures << " check(s) failed\n";
        return 1;
    }
    std::cout << "\nall checks passed\n";
    return 0;
}
