// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterTests
 * \brief Bounds on the long-wave flux coefficient as the free surface levels.
 */
#include <config.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <dumux/freeflow/shallowwater/longwave/regularization.hh>

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

using namespace Dumux::LongWave;

//! Wide-channel conveyance, independent of the model's parameter-reading version
double sheetConveyance(const double h)
{ using std::pow; return pow(h, 5.0/3.0); }

//! The flux coefficient D of q = -D grad H, with an explicit floor
double diffusivityAt(const double conveyance, const double manningN,
                     const double normGradH, const double threshold, const bool useLinear)
{ return conveyance/manningN*regularInvSqrt(normGradH, threshold, useLinear); }

/*!
 * \brief The floor that bounds D by maxDiffusivity, combined with a direct floor.
 *
 * Mirrors `gradHThreshold` without its parameter reads, which resolve once per process and
 * so cannot be swept.
 */
double thresholdFor(const double conveyance, const double manningN,
                    const double gradHEps, const double maxD, const bool useLinear)
{
    using std::max;
    return max(gradHEps, diffusivityLimitedThreshold(
        conveyance, manningN, maxD, regularInvSqrtPeak(useLinear)
    ));
}

const std::vector<double> depths{0.001, 0.01, 0.1, 1.0, 3.0, 10.0};
const std::vector<double> roughnesses{0.01, 0.025, 0.05, 0.1, 0.25};
const std::vector<double> maxDiffusivities{1e2, 1e3, 1e4, 1e5, 1e6};

void testContinuity()
{
    std::cout << "\n-- regularInvSqrt is exact above the threshold and C1 across it\n";
    for (const bool useLinear : {false, true})
    {
        const std::string tag = useLinear ? " (linear)" : " (quadratic)";
        for (const double eps : {1e-10, 1e-8, 1e-6, 1e-4})
        {
            using std::sqrt;
            for (const double factor : {1.0001, 1.1, 10.0, 1e4})
            {
                const auto x = factor*eps;
                check("exact above threshold" + tag, regularInvSqrt(x, eps, useLinear), 1.0/sqrt(x), 1e-14);
            }

            check("value continuous at threshold" + tag,
                  regularInvSqrt(eps, eps, useLinear), 1.0/sqrt(eps), 1e-14);

            // one-sided slopes at the threshold, differenced on each branch
            const auto d = 1e-6*eps;
            const auto below = (regularInvSqrt(eps, eps, useLinear)
                                - regularInvSqrt(eps - d, eps, useLinear))/d;
            const auto above = (regularInvSqrt(eps + d, eps, useLinear)
                                - regularInvSqrt(eps, eps, useLinear))/d;
            check("derivative continuous at threshold" + tag, below, above, 1e-5);
            check("derivative matches -0.5*x^(-3/2)" + tag, below, -0.5*std::pow(eps, -1.5), 1e-5);
        }
    }
}

void testPeakFactors()
{
    std::cout << "\n-- peak factors, which the diffusivity bound is derived from\n";
    checkTrue("quadratic peak is 1.25", regularInvSqrtPeak(false) == 1.25);
    checkTrue("linear peak is 1.5", regularInvSqrtPeak(true) == 1.5);

    for (const bool useLinear : {false, true})
        for (const double eps : {1e-10, 1e-8, 1e-6})
            check(std::string("peak attained at zero") + (useLinear ? " (linear)" : " (quadratic)"),
                  regularInvSqrt(0.0, eps, useLinear),
                  regularInvSqrtPeak(useLinear)/std::sqrt(eps), 1e-14);
}

void testMonotonicity()
{
    std::cout << "\n-- the flux stays monotone through the threshold, so the Jacobian cannot flip\n";
    for (const bool useLinear : {false, true})
    {
        const std::string tag = useLinear ? " (linear)" : " (quadratic)";
        const double eps = 1e-6;
        double previous = -1.0;
        bool monotone = true;
        for (int i = 0; i <= 400; ++i)
        {
            const auto x = 1e-12*std::pow(1e10, i/400.0);
            const auto flux = x*regularInvSqrt(x, eps, useLinear);
            if (flux < previous)
                monotone = false;
            previous = flux;
        }
        checkTrue("flux x/sqrt(x) non-decreasing across the threshold" + tag, monotone);
    }

    // on the quadratic branch dq/dx = 3a x^2 + c, minimal at the threshold
    const double eps = 1e-6;
    const auto a = -0.25*std::pow(eps, -2.5);
    const auto c = 1.25/std::sqrt(eps);
    check("dq/dx at the threshold is 0.5/sqrt(eps)", 3.0*a*eps*eps + c, 0.5/std::sqrt(eps), 1e-12);
}

void testCapIsRespected()
{
    std::cout << "\n-- D never exceeds MaxDiffusivity, over decades of depth, roughness and bound\n";
    bool respected = true, everBinding = false;
    double worstOvershoot = 0.0;
    for (const bool useLinear : {false, true})
        for (const auto h : depths)
            for (const auto n : roughnesses)
                for (const auto maxD : maxDiffusivities)
                {
                    const auto k = sheetConveyance(h);
                    const auto eps = thresholdFor(k, n, 0.0, maxD, useLinear);
                    for (int i = 0; i <= 200; ++i)
                    {
                        const auto x = 1e-16*std::pow(1e14, i/200.0);
                        const auto d = diffusivityAt(k, n, x, eps, useLinear);
                        if (d > maxD*(1.0 + 1e-12))
                        {
                            respected = false;
                            worstOvershoot = std::max(worstOvershoot, d/maxD - 1.0);
                        }
                        if (d > 0.99*maxD)
                            everBinding = true;
                    }
                }

    checkTrue("D <= MaxDiffusivity everywhere", respected);
    if (!respected)
        std::cerr << "      worst overshoot: " << worstOvershoot << "\n";
    checkTrue("the sweep actually reaches the bound somewhere", everBinding);
}

void testCapIsTight()
{
    std::cout << "\n-- and the bound is attained, so it does not over-restrict\n";
    for (const bool useLinear : {false, true})
        for (const auto h : depths)
            for (const auto maxD : maxDiffusivities)
            {
                const auto n = 0.03;
                const auto k = sheetConveyance(h);
                const auto eps = thresholdFor(k, n, 0.0, maxD, useLinear);
                check(std::string("sup D = MaxDiffusivity")
                          + (useLinear ? " (linear)" : " (quadratic)"),
                      diffusivityAt(k, n, 0.0, eps, useLinear), maxD, 1e-12);
            }
}

void testInertWhenNotBinding()
{
    std::cout << "\n-- an infinite bound changes nothing, and a loose one leaves sheet flow alone\n";
    const auto inf = std::numeric_limits<double>::infinity();
    for (const bool useLinear : {false, true})
        for (const auto h : depths)
            for (const auto n : roughnesses)
            {
                const auto k = sheetConveyance(h);
                checkTrue("infinite bound implies no floor of its own",
                          diffusivityLimitedThreshold(k, n, inf, regularInvSqrtPeak(useLinear)) == 0.0);
                checkTrue("infinite bound leaves GradHEpsilon as the threshold",
                          thresholdFor(k, n, 1e-8, inf, useLinear) == 1e-8);
            }

    // 1 cm of sheet flow tops out near 60 m^2/s at GradHEpsilon = 1e-8, so a bound of
    // 1e5 m^2/s must not touch it while still binding on a lake
    const auto shallow = sheetConveyance(0.01);
    const auto deep = sheetConveyance(10.0);
    checkTrue("a lake-sized bound is inert for sheet flow",
              thresholdFor(shallow, 0.1, 1e-8, 1e5, false) == 1e-8);
    checkTrue("the same bound binds on a lake",
              thresholdFor(deep, 0.03, 1e-8, 1e5, false) > 1e-8);
}

void testDepthAwareness()
{
    std::cout << "\n-- the reason the parameter exists: one slope floor is not depth-aware\n";

    // a single GradHEpsilon lets D run away with the conveyance ...
    const auto eps = 1e-8;
    const auto n = 0.03;
    const auto shallowD = diffusivityAt(sheetConveyance(0.01), n, 0.0, eps, false);
    const auto deepD = diffusivityAt(sheetConveyance(10.0), n, 0.0, eps, false);
    checkTrue("GradHEpsilon alone spans >5 decades of D between sheet flow and a lake",
              deepD/shallowD > 1e5);

    // ... while the derived floor tracks it, holding D at the bound for every depth
    const auto maxD = 1e5;
    for (const auto h : depths)
    {
        const auto k = sheetConveyance(h);
        check("D at rest is the bound regardless of depth",
              diffusivityAt(k, n, 0.0, thresholdFor(k, n, 0.0, maxD, false), false), maxD, 1e-12);
    }

    // the implied floor scales as the square of the conveyance
    const auto k1 = sheetConveyance(1.0), k2 = sheetConveyance(4.0);
    const auto t1 = diffusivityLimitedThreshold(k1, n, maxD, 1.25);
    const auto t2 = diffusivityLimitedThreshold(k2, n, maxD, 1.25);
    check("implied floor scales as conveyance^2", t2/t1, (k2/k1)*(k2/k1), 1e-12);
}

void testAdditiveShiftForm()
{
    std::cout << "\n-- the coupling interfaces use 1/sqrt(|gradH| + eps), whose peak is 1\n";
    const auto n = 0.03;
    for (const auto h : depths)
        for (const auto maxD : maxDiffusivities)
        {
            const auto k = sheetConveyance(h);
            const auto floor = diffusivityLimitedThreshold(k, n, maxD, 1.0);
            // worst case for that form is |gradH| = 0, where normH is exactly the floor
            check("sup D = MaxDiffusivity for the additive shift",
                  k/(n*std::sqrt(std::max(0.0 + 0.0, floor))), maxD, 1e-12);
        }
}

} // end anonymous namespace

int main()
{
    std::cout << "Long-wave flux coefficient bounds\n";

    testContinuity();
    testPeakFactors();
    testMonotonicity();
    testCapIsRespected();
    testCapIsTight();
    testInertWhenNotBinding();
    testDepthAwareness();
    testAdditiveShiftForm();

    if (failures > 0)
    {
        std::cerr << "\n" << failures << " check(s) failed\n";
        return 1;
    }

    std::cout << "\nAll checks passed\n";
    return 0;
}
