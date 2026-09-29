// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterTests
 * \brief Canopy interception and per-dof Green-Ampt infiltration.
 */
#include <config.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <dumux/freeflow/shallowwater/surfacerunoff/greenampt.hh>
#include <dumux/freeflow/shallowwater/surfacerunoff/infiltration.hh>

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

} // end anonymous namespace

int main()
{
    using namespace Dumux::SurfaceRunoff;
    using Soil = GreenAmptSoil<double>;
    using Losses = SurfaceLosses<double>;

    // Goodwin Creek soils and canopies: Falaya under forest, Memphis under pasture
    const Soil falaya(0.307e-2/3600.0, 0.14, 0.29, 1e6);
    const Soil memphis(0.432e-2/3600.0, 0.22, 0.29, 1e6);
    const auto forest = 3.0e-3, pasture = 1.5e-3;
    const auto rain = 40.0e-3/3600.0; // 40 mm/h, the middle of the benchmark storm
    const auto dt = 60.0;

    std::cout << "-- the canopy fills before anything reaches the ground\n";
    {
        // a soil that takes nothing isolates the canopy
        Losses losses(1, Soil(0.0, 0.0, 0.0, 0.0), forest);
        double intercepted = 0.0, throughfall = 0.0;
        const int steps = 20;
        for (int i = 0; i < steps; ++i)
        {
            losses.update(0, rain, 0.0, dt, 0.5);
            intercepted += losses.interceptionRate(0)*dt;
            throughfall += (rain - losses.lossRate(0))*dt;
            losses.commit(dt);
        }
        check("the canopy holds exactly its capacity", intercepted, forest, 1e-12);
        check("stored equals intercepted", losses.interceptionStored()[0], forest, 1e-12);
        check("the rest gets through", throughfall, rain*steps*dt - forest, 1e-12);
        // 3 mm at 40 mm/h is 4.5 min, so nothing gets through in the first four steps
        checkTrue("nothing gets through while the canopy is filling",
                  rain*4*dt <= forest);
    }

    std::cout << "\n-- interception delays infiltration rather than replacing it\n";
    {
        Losses bare(1, falaya, 0.0);
        Losses canopied(1, falaya, forest);
        const auto pond = 0.0;
        bare.update(0, rain, pond, dt, 0.5);
        canopied.update(0, rain, pond, dt, 0.5);
        checkTrue("a dry canopy takes all the rain in the first step",
                  canopied.infiltrationRate(0) < bare.infiltrationRate(0));
        // once the canopy is full the two soils see the same supply
        for (int i = 0; i < 100; ++i)
        {
            canopied.update(0, rain, pond, dt, 0.5);
            canopied.commit(dt);
            bare.update(0, rain, pond, dt, 0.5);
            bare.commit(dt);
        }
        canopied.update(0, rain, pond, dt, 0.5);
        check("the canopy is full and passes everything on", canopied.interceptionRate(0), 0.0, 1e-30);
        checkTrue("the soil under it has taken less, being 3 mm behind",
                  canopied.wettedDepth()[0] < bare.wettedDepth()[0]);
    }

    std::cout << "\n-- soils are independent per degree of freedom\n";
    {
        Losses losses(std::vector<Soil>{falaya, memphis}, std::vector<double>{forest, pasture});
        check("size", double(losses.size()), 2.0, 1e-30);
        for (int i = 0; i < 400; ++i)
        {
            for (std::size_t dof = 0; dof < 2; ++dof)
                losses.update(dof, rain, 0.0, dt, 0.5);
            losses.commit(dt);
        }
        checkTrue("the more conductive soil under the thinner canopy is wetter",
                  losses.wettedDepth()[1] > losses.wettedDepth()[0]);
        check("each canopy holds its own capacity", losses.interceptionStored()[0], forest, 1e-12);
        check("and the other holds its own", losses.interceptionStored()[1], pasture, 1e-12);
    }

    std::cout << "\n-- the stores move by the step that was actually taken\n";
    {
        Losses losses(1, falaya, 0.0);
        losses.update(0, rain, 0.0, dt, 0.5);
        const auto rate = losses.infiltrationRate(0);
        losses.commit(0.25*dt); // the solver retried at a quarter of the step
        check("wetted depth follows the committed step", losses.wettedDepth()[0], rate*0.25*dt, 1e-12);
    }

    std::cout << "\n-- losses never exceed what is available\n";
    {
        Losses losses(1, memphis, forest);
        const auto pond = 5.0e-3;
        const auto drawdown = 0.5;
        bool bounded = true;
        for (int i = 0; i < 200; ++i)
        {
            losses.update(0, rain, pond, dt, drawdown);
            bounded = bounded && losses.lossRate(0) <= rain + drawdown*pond/dt + 1e-18;
            losses.commit(dt);
        }
        checkTrue("bounded by rainfall plus the pond offered", bounded);
    }

    std::cout << "\n-- a uniform field reproduces the bare soil it was built from\n";
    {
        Losses losses(1, falaya, 0.0);
        double fieldDepth = 0.0, soilDepth = 0.0;
        for (int i = 0; i < 50; ++i)
        {
            losses.update(0, rain, 0.0, dt, 0.0);
            losses.commit(dt);
            soilDepth += falaya.infiltration(soilDepth, rain, dt);
            fieldDepth = losses.wettedDepth()[0];
        }
        check("wetted depth", fieldDepth, soilDepth, 1e-12);
    }

    std::cout << "\n-- an initially wet soil takes less\n";
    {
        Losses dry(1, falaya, 0.0), wet(1, falaya, 0.0);
        wet.setInitialWettedDepth(0.05);
        dry.update(0, rain, 0.0, dt, 0.0);
        wet.update(0, rain, 0.0, dt, 0.0);
        checkTrue("the wetter column has the lower capacity",
                  wet.infiltrationRate(0) < dry.infiltrationRate(0));
    }

    std::cout << "\n-- initial wetting can vary per degree of freedom\n";
    {
        const Soil shallow(1.0e-5, 0.1, 0.5, 0.2);
        Losses losses(std::vector<Soil>{falaya, shallow}, std::vector<double>{0.0, 0.0});
        losses.setInitialWettedDepth(std::vector<double>{0.02, 1.0});
        check("the first column keeps its prescribed depth", losses.wettedDepth()[0], 0.02, 1e-12);
        check("the second column is capped by its store", losses.wettedDepth()[1], shallow.capacity(), 1e-12);

        Losses dry(2, falaya, 0.0), wet(2, falaya, 0.0);
        wet.setInitialWettedDepth(std::vector<double>{0.02, 0.08});
        for (std::size_t dof = 0; dof < wet.size(); ++dof)
        {
            dry.update(dof, rain, 0.0, dt, 0.0);
            wet.update(dof, rain, 0.0, dt, 0.0);
        }
        checkTrue("the wetter prescribed column has the lower capacity",
                  wet.infiltrationRate(1) < wet.infiltrationRate(0));
        checkTrue("both prescribed columns start below the dry reference capacity",
                  wet.infiltrationRate(0) < dry.infiltrationRate(0)
                  && wet.infiltrationRate(1) < dry.infiltrationRate(1));
    }

    if (failures > 0)
    {
        std::cerr << "\n" << failures << " check(s) failed\n";
        return 1;
    }
    std::cout << "\nall checks passed\n";
    return 0;
}
