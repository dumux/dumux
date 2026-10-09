// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkModels
 * \brief Test of the pore and throat geometry of a regular triangulation of a sphere packing
 */
#include <config.h>

#include <array>
#include <cmath>
#include <iostream>
#include <numbers>
#include <random>
#include <string>

#include <dune/common/exceptions.hh>

#include <dumux/porenetwork/extraction/spherepackinggeometry.hh>

namespace {

using namespace Dumux::PoreNetwork::SpherePacking;
using P = Point<double>;
constexpr double pi = std::numbers::pi;

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

std::array<P, 4> regularTetrahedron(double a)
{
    const double s = a/std::sqrt(8.0);
    return {P{s, s, s}, P{s, -s, -s}, P{-s, s, -s}, P{-s, -s, s}};
}

std::array<P, 3> equilateralTriangle(double a)
{
    const double h = a*std::sqrt(3.0)/2.0;
    return {P{0.0, 0.0, 0.0}, P{a, 0.0, 0.0}, P{0.5*a, h, 0.0}};
}

bool inTetrahedron(const P& q, const std::array<P, 4>& x)
{
    const double v = tetrahedronVolume(x);
    double sum = 0.0;
    for (int i = 0; i < 4; ++i)
    {
        auto y = x;
        y[i] = q;
        sum += tetrahedronVolume(y);
    }
    return sum <= v*(1.0 + 1e-12);
}

bool inAnySphere(const P& q, const auto& x, const auto& r)
{
    for (std::size_t i = 0; i < x.size(); ++i)
        if ((q - x[i]).two_norm2() < r[i]*r[i])
            return true;
    return false;
}

void testSolidAngle()
{
    const P o{0.0, 0.0, 0.0};
    check(solidAngle(o, P{1.0, 0.0, 0.0}, P{0.0, 1.0, 0.0}, P{0.0, 0.0, 1.0}), pi/2.0, 1e-14, "solid angle of a cube corner");

    const auto x = regularTetrahedron(1.0);
    check(solidAngle(x[0], x[1], x[2], x[3]), std::acos(23.0/27.0), 1e-14, "solid angle of a regular tetrahedron corner");

    // apex slightly above the centroid of a large triangle sees almost a half space
    const auto t = equilateralTriangle(1.0);
    const P centroid = (t[0] + t[1] + t[2])/3.0;
    const double height = 1e-3;
    const double omega = solidAngle(centroid + P{0.0, 0.0, height}, t[0], t[1], t[2]);
    if (!(omega > pi && omega < 2.0*pi))
    {
        std::cout << "FAILED solid angle above pi: " << omega << "\n";
        ++numFailures;
    }
    // the three triangles (apex, edge) and the base close the tetrahedron: their solid angles seen from
    // a point inside sum to 4 pi
    const auto y = regularTetrahedron(1.3);
    const P q{0.05, -0.02, 0.07};
    double sum = 0.0;
    for (int i = 0; i < 4; ++i)
        sum += solidAngle(q, y[(i+1)%4], y[(i+2)%4], y[(i+3)%4]);
    check(sum, 4.0*pi, 1e-13, "solid angles around an interior point");
}

void testRegularConfigurations()
{
    const double a = 2.0, r = 1.0;
    const auto x = regularTetrahedron(a);
    const std::array<double, 4> radii{r, r, r, r};

    const auto c = powerCenter(x, radii);
    check(c.two_norm() + 1.0, 1.0, 1e-14, "power centre of the regular tetrahedron at its centroid");

    const auto sphere = inscribedSphere(x, radii);
    check(sphere.radius, a*std::sqrt(6.0)/4.0 - r, 1e-14, "inscribed sphere of the regular tetrahedron");

    const double voidVolume = a*a*a/(6.0*std::sqrt(2.0)) - 4.0*std::acos(23.0/27.0)*r*r*r/3.0;
    check(tetrahedronVoidVolume(x, radii), voidVolume, 1e-14, "void volume of the regular tetrahedron");

    const auto t = equilateralTriangle(a);
    const std::array<double, 3> r3{r, r, r};
    const auto circle = inscribedCircle(t, r3);
    check(circle.radius, a/std::sqrt(3.0) - r, 1e-14, "inscribed circle of the equilateral triangle");
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
        const double power0 = (c - x[0]).two_norm2() - r[0]*r[0];
        for (int i = 1; i < 4; ++i)
            check((c - x[i]).two_norm2() - r[i]*r[i], power0, 1e-11, "equal power of the power centre");

        const auto s = inscribedSphere(x, r);
        if (s.radius > 0.0)
            for (int i = 0; i < 4; ++i)
                check((s.center - x[i]).two_norm(), r[i] + s.radius, 1e-12, "tangency of the inscribed sphere");

        const std::array<P, 3> t{x[0], x[1], x[2]};
        const std::array<double, 3> rt{r[0], r[1], r[2]};
        const auto circle = inscribedCircle(t, rt);
        if (circle.radius > 0.0)
        {
            for (int i = 0; i < 3; ++i)
                check((circle.center - t[i]).two_norm(), rt[i] + circle.radius, 1e-12, "tangency of the inscribed circle");
            const auto normal = Dumux::PoreNetwork::SpherePacking::Detail::cross(t[1] - t[0], t[2] - t[0]);
            check(1.0 + (circle.center - t[0])*normal/normal.two_norm(), 1.0, 1e-13, "inscribed circle in the facet plane");
        }
    }
}

// configuration in which no sphere reaches beyond the opposite face, so the sector formulas are exact
std::pair<std::array<P, 4>, std::array<double, 4>> randomTetrahedron()
{
    std::array<P, 4> x = regularTetrahedron(1.0);
    x[0] += P{0.05, -0.08, 0.03};
    x[2] += P{-0.06, 0.02, 0.07};
    return {x, {0.42, 0.36, 0.47, 0.39}};
}

void testMonteCarlo()
{
    std::mt19937 gen(7);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    const int n = 4'000'000;

    {
        const auto [x, r] = randomTetrahedron();
        P lower(1e9), upper(-1e9);
        for (const auto& v : x)
            for (int k = 0; k < 3; ++k)
            {
                lower[k] = std::min(lower[k], v[k]);
                upper[k] = std::max(upper[k], v[k]);
            }
        const auto extent = upper - lower;
        int hits = 0;
        for (int s = 0; s < n; ++s)
        {
            const P q{lower[0] + extent[0]*unit(gen), lower[1] + extent[1]*unit(gen), lower[2] + extent[2]*unit(gen)};
            if (inTetrahedron(q, x) && !inAnySphere(q, x, r))
                ++hits;
        }
        const double estimate = double(hits)/n*extent[0]*extent[1]*extent[2];
        check(tetrahedronVoidVolume(x, r), estimate, 5e-3, "void volume vs Monte Carlo");
    }

    {
        // circle 0 crosses the opposite edge, so the circular segment beyond it is added back
        const std::array<P, 3> t{P{0.0, 0.0, 0.0}, P{-1.0, 0.4, 0.0}, P{1.2, 0.4, 0.0}};
        const std::array<double, 3> r{0.5, 0.3, 0.35};
        if (circularSegmentBeyondEdge(t[0], r[0], t[1], t[2]) <= 0.0)
        {
            std::cout << "FAILED circle does not cross the opposite edge\n";
            ++numFailures;
        }
        int hits = 0, inside = 0;
        for (int s = 0; s < n; ++s)
        {
            const P q{-1.0 + 2.2*unit(gen), 0.4*unit(gen), 0.0};
            const double u = q[1]/0.4;
            if (q[0] < -u || q[0] > 1.2*u)
                continue;
            ++inside;
            if (!inAnySphere(q, t, r))
                ++hits;
        }
        const double area = triangleArea(t[0], t[1], t[2]);
        check(facetFluidArea(t, r), area*hits/inside, 5e-3, "fluid area with a circle crossing an edge vs Monte Carlo");
    }

    {
        // throat region of a facet between two power centres on opposite sides
        const std::array<P, 3> t = equilateralTriangle(1.0);
        const std::array<double, 3> r{0.42, 0.38, 0.45};
        const P centroid = (t[0] + t[1] + t[2])/3.0;
        const P p1 = centroid + P{0.03, -0.02, 0.35};
        const P p2 = centroid + P{-0.01, 0.04, -0.28};
        const auto region = throatRegion(t, r, p1, p2);

        const std::array<P, 4> upper{t[0], t[1], t[2], p1};
        const std::array<P, 4> lowerTet{t[0], t[1], t[2], p2};
        check(region.volume, tetrahedronVolume(upper) + tetrahedronVolume(lowerTet), 1e-13, "bipyramid volume");

        int hits = 0;
        for (int s = 0; s < n; ++s)
        {
            const P q{1.0*unit(gen), 0.9*unit(gen), -0.3 + 0.7*unit(gen)};
            if ((inTetrahedron(q, upper) || inTetrahedron(q, lowerTet)) && !inAnySphere(q, t, r))
                ++hits;
        }
        check(region.voidVolume, double(hits)/n*1.0*0.9*0.7, 5e-3, "throat void volume vs Monte Carlo");

        std::normal_distribution<double> normal;
        double surface = 0.0;
        for (int i = 0; i < 3; ++i)
        {
            int surfaceHits = 0;
            for (int s = 0; s < n/4; ++s)
            {
                P dir{normal(gen), normal(gen), normal(gen)};
                const P q = t[i] + dir*(r[i]/dir.two_norm());
                if (inTetrahedron(q, upper) || inTetrahedron(q, lowerTet))
                    ++surfaceHits;
            }
            surface += double(surfaceHits)/(n/4)*4.0*pi*r[i]*r[i];
        }
        check(region.solidSurface, surface, 5e-3, "throat solid surface vs Monte Carlo");
    }
}

} // end anonymous namespace

int main()
{
    testSolidAngle();
    testRegularConfigurations();
    testRandomTangency();
    testMonteCarlo();

    if (numFailures > 0)
        DUNE_THROW(Dune::Exception, numFailures << " checks failed");
    std::cout << "all checks passed\n";
    return 0;
}
