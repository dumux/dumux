// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterTests
 * \brief Stage-discharge relations for channel outlets.
 */
#include <config.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <dumux/freeflow/shallowwater/surfacerunoff/outletcontrol.hh>

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

using namespace Dumux::SurfaceRunoff;

const std::vector<double> coefficients{1.4, 1.7, 1.84, 2.0};
const std::vector<double> crests{-6.07, -1.0, 0.0, 2.598, 124.0};

void testWeirLaw()
{
    std::cout << "\n-- the weir law is C*head^(3/2), and dry below the crest\n";
    for (const auto c : coefficients)
        for (const auto crest : crests)
        {
            // the tolerance is not machine precision because the head cannot be: forming
            // it as (crest + head) - crest cancels, and at a crest of 124 m a 0.1 mm head
            // survives to only ~1e-12 relative. That is a property of absolute elevations,
            // not of the law, and it is why the header warns about it.
            for (const auto head : {1e-4, 1e-3, 0.01, 0.1, 1.0, 3.0})
                check("q = C*head^1.5", weirDischargePerLength(crest + head, crest, c),
                      c*std::pow((crest + head) - crest, 1.5), 1e-12);

            checkTrue("dry at the crest", weirDischargePerLength(crest, crest, c) == 0.0);
            for (const auto below : {1e-6, 0.01, 1.0, 5.0})
                checkTrue("dry below the crest",
                          weirDischargePerLength(crest - below, crest, c) == 0.0);
        }
}

void testInverse()
{
    std::cout << "\n-- the inverse recovers the head, over decades of discharge\n";
    for (const auto c : coefficients)
        for (int i = 0; i <= 30; ++i)
        {
            const auto head = 1e-5*std::pow(1e6, i/30.0);
            const auto q = weirDischargePerLength(head, 0.0, c);
            check("head -> q -> head", weirHeadForDischargePerLength(q, c), head, 1e-12);
        }

    std::cout << "\n-- and a non-positive discharge means no head, not a domain error\n";
    for (const auto c : coefficients)
        for (const auto q : {0.0, -1e-12, -1.0})
            checkTrue("no head below zero discharge",
                      weirHeadForDischargePerLength(q, c) == 0.0);
}

void testSmoothnessAtCrest()
{
    std::cout << "\n-- C1 at the crest, so Newton sees no kink where the weir starts to flow\n";
    for (const auto c : coefficients)
    {
        const auto crest = 2.598;
        const auto d = 1e-7;

        // one-sided slopes: dry below, 1.5*C*sqrt(head) above, both vanishing at the crest
        const auto below = (weirDischargePerLength(crest, crest, c)
                            - weirDischargePerLength(crest - d, crest, c))/d;
        const auto above = (weirDischargePerLength(crest + d, crest, c)
                            - weirDischargePerLength(crest, crest, c))/d;
        checkTrue("slope vanishes below the crest", below == 0.0);
        checkTrue("slope vanishes at the crest from above", above < 1e-3);

        // and away from the crest the slope is the analytic one
        for (const auto head : {0.01, 0.5, 2.0})
        {
            const auto numeric = (weirDischargePerLength(crest + head + d, crest, c)
                                  - weirDischargePerLength(crest + head - d, crest, c))/(2*d);
            check("dq/dH = 1.5*C*sqrt(head)", numeric, 1.5*c*std::sqrt(head), 1e-6);
        }
    }
}

void testMonotonicity()
{
    std::cout << "\n-- discharge rises monotonically with stage\n";
    for (const auto c : coefficients)
    {
        bool monotone = true;
        double previous = -1.0;
        for (int i = 0; i <= 500; ++i)
        {
            const auto level = -1.0 + 4.0*i/500.0;
            const auto q = weirDischargePerLength(level, 0.0, c);
            if (q < previous)
                monotone = false;
            previous = q;
        }
        checkTrue("q non-decreasing in stage across the crest", monotone);
    }
}

void testNames()
{
    std::cout << "\n-- names round-trip, and an unknown one is rejected rather than ignored\n";
    for (const auto control : {OutletControl::closed, OutletControl::freeOutfall,
                               OutletControl::weir, OutletControl::prescribedDischarge,
                               OutletControl::fixedStage})
        checkTrue("name round-trips",
                  outletControlFromName(outletControlName(control)) == control);

    bool threw = false;
    try { outletControlFromName("dam"); }
    catch (const Dumux::ParameterException&) { threw = true; }
    checkTrue("an unknown control throws", threw);
}

/*!
 * \brief The property the lake benchmark exists to check, in closed form.
 *
 * A lake under steady inflow settles where the weir passes exactly what arrives. Solving
 * for that level needs no simulation, so the benchmark has a reference to be scored against
 * rather than merely a plausible-looking hydrograph.
 */
void testSteadyStatePrediction()
{
    std::cout << "\n-- steady lake level under a steady inflow, from the weir law alone\n";
    const auto crestLength = 20.0;
    for (const auto c : coefficients)
        for (const auto crest : crests)
            for (const auto inflow : {5.05e-3, 0.1, 1.0, 25.0})
            {
                const auto level = crest + weirHeadForDischargePerLength(inflow/crestLength, c);
                const auto outflow = crestLength*weirDischargePerLength(level, crest, c);
                check("weir passes exactly the inflow at the predicted level", outflow, inflow, 1e-10);
            }
}

} // end anonymous namespace

int main()
{
    std::cout << "Outlet stage-discharge relations\n";

    testWeirLaw();
    testInverse();
    testSmoothnessAtCrest();
    testMonotonicity();
    testNames();
    testSteadyStatePrediction();

    if (failures > 0)
    {
        std::cerr << "\n" << failures << " check(s) failed\n";
        return 1;
    }

    std::cout << "\nAll checks passed\n";
    return 0;
}
