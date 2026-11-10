// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkModels
 * \brief Test of the pore and throat geometry of a regular triangulation of a sphere packing: exact values for
 *        regular configurations, tangency of the inscribed spheres and circles, and Monte Carlo estimates of void
 *        volumes, fluid areas and throat regions
 */
#include <config.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <numbers>
#include <random>
#include <string>

#include <dune/common/exceptions.hh>

#include <dumux/common/math.hh>
#include <dumux/porenetwork/extraction/spherepackinggeometry.hh>

namespace {

using namespace Dumux::PoreNetwork::SpherePacking;
using P = Point<double>;
constexpr double pi = std::numbers::pi;
constexpr int numSamples = 4'000'000;
// a correct formula fails a Monte Carlo check with probability 6e-5
constexpr double numStandardErrors = 4.0;

int numFailures = 0;

void check(double value, double reference, double tolerance, const std::string& name)
{
    const double error = std::abs(value - reference)/std::max(std::abs(reference), 1e-300);
    if (!(error <= tolerance))
    {
        std::cout << "FAILED " << name << ": " << value << " vs " << reference << " (relative error " << error << ")\n";
        ++numFailures;
    }
}

void checkTrue(bool condition, const std::string& name)
{
    if (!condition)
    {
        std::cout << "FAILED " << name << "\n";
        ++numFailures;
    }
}

std::array<P, 4> regularTetrahedron(double a)
{
    const double s = a/std::sqrt(8.0);
    return {P{s, s, s}, P{s, -s, -s}, P{-s, s, -s}, P{-s, -s, s}};
}

std::array<P, 3> equilateralTriangle(double a)
{ return {P{0.0, 0.0, 0.0}, P{a, 0.0, 0.0}, P{0.5*a, a*std::sqrt(3.0)/2.0, 0.0}}; }

bool inTetrahedron(const P& q, const std::array<P, 4>& x)
{
    double sum = 0.0;
    for (int i = 0; i < 4; ++i)
    {
        auto y = x;
        y[i] = q;
        sum += tetrahedronVolume(y);
    }
    return sum <= tetrahedronVolume(x)*(1.0 + 1e-12);
}

bool inAnySphere(const P& q, const auto& x, const auto& r)
{
    for (std::size_t i = 0; i < x.size(); ++i)
        if ((q - x[i]).two_norm2() < r[i]*r[i])
            return true;
    return false;
}

//! integral of f over the box [lower, upper] (over a rectangle if one extent is zero), with its standard error
std::array<double, 2> integral(const P& lower, const P& upper, const std::function<double(const P&)>& f, std::mt19937& gen)
{
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    double measure = 1.0, sum = 0.0, sumSquares = 0.0;
    for (int k = 0; k < 3; ++k)
        measure *= upper[k] > lower[k] ? upper[k] - lower[k] : 1.0;
    for (int s = 0; s < numSamples; ++s)
    {
        P q;
        for (int k = 0; k < 3; ++k)
            q[k] = lower[k] + (upper[k] - lower[k])*unit(gen);
        const double value = f(q);
        sum += value;
        sumSquares += value*value;
    }
    const double mean = sum/numSamples;
    return {mean*measure, std::sqrt((sumSquares/numSamples - mean*mean)/numSamples)*measure};
}

//! integral of f over the surface of a sphere, with its standard error
std::array<double, 2> surfaceIntegral(const P& center, double radius, const std::function<double(const P&)>& f, std::mt19937& gen)
{
    std::normal_distribution<double> normal;
    double sum = 0.0, sumSquares = 0.0;
    for (int s = 0; s < numSamples; ++s)
    {
        P direction{normal(gen), normal(gen), normal(gen)};
        const double value = f(center + direction*(radius/direction.two_norm()));
        sum += value;
        sumSquares += value*value;
    }
    const double mean = sum/numSamples, area = 4.0*pi*radius*radius;
    return {mean*area, std::sqrt((sumSquares/numSamples - mean*mean)/numSamples)*area};
}

void checkEstimate(double value, const std::array<double, 2>& estimate, const std::string& name)
{ check(value, estimate[0], numStandardErrors*estimate[1]/std::abs(estimate[0]), name + " vs Monte Carlo"); }

void testRegularConfigurations()
{
    const P o(0.0);
    check(solidAngle(o, P{1.0, 0.0, 0.0}, P{0.0, 1.0, 0.0}, P{0.0, 0.0, 1.0}), pi/2.0, 1e-14, "solid angle of a cube corner");
    const auto t0 = equilateralTriangle(1.0);
    const P above = (t0[0] + t0[1] + t0[2])/3.0 + P{0.0, 0.0, 1e-3};
    const double halfSpace = solidAngle(above, t0[0], t0[1], t0[2]);
    checkTrue(halfSpace > pi && halfSpace < 2.0*pi, "solid angle above pi");

    const double a = 2.0, r = 1.0;
    const auto x = regularTetrahedron(a);
    const std::array<double, 4> radii{r, r, r, r};
    check(1.0 + powerCenter(x, radii).two_norm(), 1.0, 1e-14, "power centre of the regular tetrahedron at its centroid");
    check(inscribedSphere(x, radii).radius, a*std::sqrt(6.0)/4.0 - r, 1e-14, "inscribed sphere of the regular tetrahedron");
    check(tetrahedronVoidVolume(x, radii), a*a*a/(6.0*std::sqrt(2.0)) - 4.0*std::acos(23.0/27.0)*r*r*r/3.0, 1e-14,
          "void volume of the regular tetrahedron");

    const auto t = equilateralTriangle(a);
    const std::array<double, 3> r3{r, r, r};
    check(inscribedCircle(t, r3).radius, a/std::sqrt(3.0) - r, 1e-14, "inscribed circle of the equilateral triangle");
    check(facetFluidArea(t, r3), std::sqrt(3.0)/4.0*a*a - 0.5*pi*r*r, 1e-14, "fluid area of the equilateral facet");
}

void testRandomTangency()
{
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> shift(-0.25, 0.25), radius(0.4, 0.6);
    for (int sample = 0; sample < 1000; ++sample)
    {
        auto x = regularTetrahedron(1.2);
        std::array<double, 4> r;
        for (int i = 0; i < 4; ++i)
        {
            for (int k = 0; k < 3; ++k)
                x[i][k] += shift(gen);
            r[i] = radius(gen);
        }

        const auto c = powerCenter(x, r);
        for (int i = 1; i < 4; ++i)
            check((c - x[i]).two_norm2() - r[i]*r[i], (c - x[0]).two_norm2() - r[0]*r[0], 1e-11, "equal power of the power centre");
        if (const auto s = inscribedSphere(x, r); s.radius > 0.0)
            for (int i = 0; i < 4; ++i)
                check((s.center - x[i]).two_norm(), r[i] + s.radius, 1e-12, "tangency of the inscribed sphere");

        const std::array<P, 3> t{x[0], x[1], x[2]};
        if (const auto circle = inscribedCircle(t, std::array<double, 3>{r[0], r[1], r[2]}); circle.radius > 0.0)
        {
            for (int i = 0; i < 3; ++i)
                check((circle.center - t[i]).two_norm(), r[i] + circle.radius, 1e-12, "tangency of the inscribed circle");
            check(1.0 + Dumux::tripleProduct(circle.center - t[0], t[1] - t[0], t[2] - t[0]), 1.0, 1e-13,
                  "inscribed circle in the facet plane");
        }
    }
}

//! the region of a throat equals the pyramids over the facet with apexes p1, p2, weighted with the signs s1, s2
void checkThroatRegion(const std::array<P, 3>& t, const std::array<double, 3>& r, const P& p1, const P& p2,
                       double s1, double s2, const std::string& name, std::mt19937& gen)
{
    const auto region = throatRegion(t, r, p1, p2);
    const std::array<P, 4> pyramid1{t[0], t[1], t[2], p1}, pyramid2{t[0], t[1], t[2], p2};
    check(region.volume, s1*tetrahedronVolume(pyramid1) + s2*tetrahedronVolume(pyramid2), 1e-13, name + ": volume");
    const auto indicator = [&](const P& q) { return s1*inTetrahedron(q, pyramid1) + s2*inTetrahedron(q, pyramid2); };

    P lower(1e9), upper(-1e9);
    for (const auto& v : {t[0], t[1], t[2], p1, p2})
        for (int k = 0; k < 3; ++k)
        {
            lower[k] = std::min(lower[k], v[k]);
            upper[k] = std::max(upper[k], v[k]);
        }
    checkEstimate(region.voidVolume, integral(lower, upper, [&](const P& q) { return inAnySphere(q, t, r) ? 0.0 : indicator(q); }, gen),
                  name + ": void volume");

    double surface = 0.0, variance = 0.0;
    for (int i = 0; i < 3; ++i)
    {
        const auto [value, error] = surfaceIntegral(t[i], r[i], indicator, gen);
        surface += value;
        variance += error*error;
    }
    checkEstimate(region.solidSurface, {surface, std::sqrt(variance)}, name + ": wetted surface");
}

void testMonteCarlo()
{
    std::mt19937 gen(7);

    // no sphere reaches beyond the opposite face, so the sector formulas are exact
    auto x = regularTetrahedron(1.0);
    x[0] += P{0.05, -0.08, 0.03};
    x[2] += P{-0.06, 0.02, 0.07};
    const std::array<double, 4> r{0.42, 0.36, 0.47, 0.39};
    checkEstimate(tetrahedronVoidVolume(x, r), integral(P(-0.5), P(0.5), [&](const P& q) {
        return inTetrahedron(q, x) && !inAnySphere(q, x, r); }, gen), "void volume of a tetrahedron");

    // circle 0 crosses the opposite edge, so the circular segment beyond it is added back
    const std::array<P, 3> crossing{P{0.0, 0.0, 0.0}, P{-1.0, 0.4, 0.0}, P{1.2, 0.4, 0.0}};
    const std::array<double, 3> rc{0.5, 0.3, 0.35};
    checkTrue(circularSegmentBeyondEdge(crossing[0], rc[0], crossing[1], crossing[2]) > 0.0, "circle crosses the opposite edge");
    checkEstimate(facetFluidArea(crossing, rc), integral(P{-1.0, 0.0, 0.0}, P{1.2, 0.4, 0.0}, [&](const P& q) {
        const double u = q[1]/0.4;
        return q[0] >= -u && q[0] <= 1.2*u && !inAnySphere(q, crossing, rc); }, gen), "fluid area with a circle crossing an edge");

    // throat regions: power centres on opposite sides of the facet, and on one side of it on a normal through a
    // point outside the facet, where the pyramid of the nearer power centre counts negatively
    const auto t = equilateralTriangle(1.0);
    const std::array<double, 3> rt{0.42, 0.38, 0.45};
    const P centroid = (t[0] + t[1] + t[2])/3.0, foot{0.55, -0.12, 0.0};
    checkThroatRegion(t, rt, centroid + P{0.03, -0.02, 0.35}, centroid + P{-0.01, 0.04, -0.28}, 1.0, 1.0,
                      "throat region, power centres on opposite sides", gen);
    checkThroatRegion(t, rt, foot + P{0.0, 0.0, 0.4}, foot + P{0.0, 0.0, 0.15}, 1.0, -1.0,
                      "throat region, power centres on one side", gen);
    check(throatRegion(t, rt, foot + P{0.0, 0.0, 0.15}, foot + P{0.0, 0.0, 0.4}).voidVolume,
          throatRegion(t, rt, foot + P{0.0, 0.0, 0.4}, foot + P{0.0, 0.0, 0.15}).voidVolume, 0.0,
          "throat region symmetric in the power centres");
}

} // end anonymous namespace

int main()
{
    testRegularConfigurations();
    testRandomTangency();
    testMonteCarlo();

    if (numFailures > 0)
        DUNE_THROW(Dune::Exception, numFailures << " checks failed");
    std::cout << "all checks passed\n";
    return 0;
}
