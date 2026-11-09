// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkModels
 * \brief Test of the pore network of a sphere packing with edge bodies (overlap lenses, rings) and solid tetrahedra
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
#include <vector>

#include <dune/common/exceptions.hh>

#include <dumux/common/initialize.hh>
#include <dumux/porenetwork/extraction/spherepackingedgebodies.hh>
#include <dumux/porenetwork/extraction/spherepackingtriangulation.hh>

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

void checkTrue(bool condition, const std::string& name)
{
    if (!condition)
    {
        std::cout << "FAILED " << name << "\n";
        ++numFailures;
    }
}

//! Monte Carlo estimate of the volume of the points of a box for which inside holds, with its standard error
std::array<double, 2> monteCarloVolume(const P& lower, const P& upper, const std::function<bool(const P&)>& inside,
                                       long samples, std::mt19937_64& rng)
{
    std::uniform_real_distribution<double> u(0.0, 1.0);
    long hits = 0;
    for (long s = 0; s < samples; ++s)
    {
        P x;
        for (int k = 0; k < 3; ++k)
            x[k] = lower[k] + u(rng)*(upper[k] - lower[k]);
        hits += inside(x);
    }
    const double boxVolume = (upper[0] - lower[0])*(upper[1] - lower[1])*(upper[2] - lower[2]);
    const double p = static_cast<double>(hits)/samples;
    return {p*boxVolume, std::sqrt(p*(1.0 - p)/samples)*boxVolume};
}

bool inTetrahedron(const P& x, const std::array<P, 4>& v)
{
    const auto sign = [](const P& a, const P& b, const P& c, const P& d) {
        return Detail::tripleProduct(b - a, c - a, d - a);
    };
    const double s = sign(v[0], v[1], v[2], v[3]);
    return sign(x, v[1], v[2], v[3])*s >= 0.0 && sign(v[0], x, v[2], v[3])*s >= 0.0
           && sign(v[0], v[1], x, v[3])*s >= 0.0 && sign(v[0], v[1], v[2], x)*s >= 0.0;
}

bool inBall(const P& x, const P& c, double r)
{ return (x - c).two_norm2() < r*r; }

//! solid torus with axis through a and b, centred at the point a + axial (b - a)/|b - a|
struct Torus
{
    P a, b;
    double axial, major, minor;

    bool contains(const P& x) const
    {
        auto e = b - a;
        e /= e.two_norm();
        const auto d = x - a;
        const double z = d*e - axial;
        auto radial = d;
        radial.axpy(-(d*e), e);
        const double rho = radial.two_norm() - major;
        return rho*rho + z*z < minor*minor;
    }

    P lower() const { const auto c = center(); return P{c[0] - major - minor, c[1] - major - minor, c[2] - major - minor}; }
    P upper() const { const auto c = center(); return P{c[0] + major + minor, c[1] + major + minor, c[2] + major + minor}; }
    P center() const { auto e = b - a; e /= e.two_norm(); auto c = a; c.axpy(axial, e); return c; }
    EdgeBody<double> body() const
    { return {2.0*pi*pi*major*minor*minor, pi*minor*minor, 4.0*pi*pi*major*minor}; }
};

void testLens(std::mt19937_64& rng)
{
    const double r1 = 0.6, r2 = 0.45, d = 0.9;
    const auto lens = sphereLens(r1, r2, d);
    const P c1{0.0, 0.0, 0.0}, c2{d, 0.0, 0.0};
    const auto [volume, error] = monteCarloVolume(P{d - r2, -r2, -r2}, P{r1, r2, r2},
        [&](const P& x) { return inBall(x, c1, r1) && inBall(x, c2, r2); }, 4000000, rng);
    check(lens.volume, volume, 5.0*error/volume, "lens volume vs Monte Carlo");
    const auto [area, areaError] = monteCarloVolume(P{d - r2, -r2, 0.0}, P{r1, r2, 1.0},
        [&](const P& x) { const P y{x[0], x[1], 0.0}; return inBall(y, c1, r1) && inBall(y, c2, r2); }, 4000000, rng);
    check(lens.area, area, 5.0*areaError/area, "lens area vs Monte Carlo");
    check(lens.radius*lens.radius, r1*r1 - lens.radicalDistance*lens.radicalDistance, 1e-14, "lens radius");
    check(lens.radius*lens.radius, r2*r2 - (d - lens.radicalDistance)*(d - lens.radicalDistance), 1e-13, "lens radius from sphere 2");

    const double r = 0.52, h = r - 0.5;
    check(sphereLens(r, r, 1.0).volume, 2.0*pi*h*h*(3.0*r - h)/3.0, 1e-13, "lens of equal spheres");
    checkTrue(sphereLens(0.4, 0.4, 1.0).volume == 0.0, "no lens of separated spheres");
}

/*!
 * Each correction is the body's quantity times the dihedral share; checked by sampling near the body only.
 * Regular tetrahedron of edge 1: spheres 0 and 1 overlap (lens on edge 01), a torus sits in the gap of
 * edge 23.
 */
void testShares(std::mt19937_64& rng)
{
    const std::array<P, 4> x{P{0.0, 0.0, 0.0}, P{1.0, 0.0, 0.0}, P{0.5, std::sqrt(3.0)/2.0, 0.0},
                             P{0.5, std::sqrt(3.0)/6.0, std::sqrt(2.0/3.0)}};
    const std::array<double, 4> r{0.51, 0.51, 0.3, 0.3};
    const Torus torus{x[2], x[3], 0.5, 0.05, 0.02};
    const auto lens = sphereLens(r[0], r[1], 1.0);
    const auto lensBody = overlapCorrection(r[0], r[1], 1.0);
    const double dihedral01 = dihedralAngle(x[0], x[1], x[2], x[3]);
    const double dihedral23 = dihedralAngle(x[2], x[3], x[0], x[1]);
    check(dihedral01, std::acos(1.0/3.0), 1e-14, "dihedral angle of the regular tetrahedron");

    // tetrahedron
    const P lensLower{1.0 - r[1], -lens.radius, -lens.radius}, lensUpper{r[0], lens.radius, lens.radius};
    const auto [lensInCell, lensError] = monteCarloVolume(lensLower, lensUpper,
        [&](const P& p) { return inTetrahedron(p, x) && inBall(p, x[0], r[0]) && inBall(p, x[1], r[1]); }, 4000000, rng);
    check(-lensBody.solidVolume*dihedral01/(2.0*pi), lensInCell, 5.0*lensError/lensInCell, "lens part in the tetrahedron");
    const auto [torusInCell, torusError] = monteCarloVolume(torus.lower(), torus.upper(),
        [&](const P& p) { return inTetrahedron(p, x) && torus.contains(p); }, 4000000, rng);
    check(torus.body().solidVolume*dihedral23/(2.0*pi), torusInCell, 5.0*torusError/torusInCell, "torus part in the tetrahedron");

    // facet (0, 1, 2) holds half the lens section, facet (0, 2, 3) half the torus section
    const auto inTriangle = [](const P& p, const P& a, const P& b, const P& c) {
        const auto n = Detail::cross(b - a, c - a);
        return Detail::cross(b - a, p - a)*n >= 0.0 && Detail::cross(c - b, p - b)*n >= 0.0
               && Detail::cross(a - c, p - c)*n >= 0.0;
    };
    const auto [lensInFacet, lensFacetError] = monteCarloVolume(P{1.0 - r[1], -lens.radius, 0.0}, P{r[0], lens.radius, 1.0},
        [&](const P& p) { const P y{p[0], p[1], 0.0}; return inTriangle(y, x[0], x[1], x[2]) && inBall(y, x[0], r[0]) && inBall(y, x[1], r[1]); },
        4000000, rng);
    check(-lensBody.halfSection, lensInFacet, 5.0*lensFacetError/lensInFacet, "half lens section in the facet");

    // throat region of facet (0, 1, 2) between apexes on both sides, and of facet (2, 3, 0)
    const auto bipyramid = [&](const P& p, const P& a, const P& b, const P& c, const P& p1, const P& p2) {
        return inTetrahedron(p, {a, b, c, p1}) || inTetrahedron(p, {a, b, c, p2});
    };
    const P apex1 = (x[0] + x[1] + x[2] + x[3])/4.0;
    P apex2{apex1[0], apex1[1], -0.7*apex1[2]};
    const double regionShare01 = dihedralAngle(x[0], x[1], apex1, apex2)/(2.0*pi);
    const auto [lensInRegion, lensRegionError] = monteCarloVolume(lensLower, lensUpper,
        [&](const P& p) { return bipyramid(p, x[0], x[1], x[2], apex1, apex2) && inBall(p, x[0], r[0]) && inBall(p, x[1], r[1]); },
        4000000, rng);
    check(-lensBody.solidVolume*regionShare01, lensInRegion, 5.0*lensRegionError/lensInRegion, "lens part in the throat region");

    // the facet (2, 3, 0) and apexes on both sides of it
    const auto n230 = Detail::cross(x[3] - x[2], x[0] - x[2]);
    const auto mid230 = (x[2] + x[3] + x[0])/3.0;
    P apex3 = mid230, apex4 = mid230;
    apex3.axpy(0.2/n230.two_norm(), n230);
    apex4.axpy(-0.3/n230.two_norm(), n230);
    const double regionShare23 = dihedralAngle(x[2], x[3], apex3, apex4)/(2.0*pi);
    const auto [torusInRegion, torusRegionError] = monteCarloVolume(torus.lower(), torus.upper(),
        [&](const P& p) { return bipyramid(p, x[2], x[3], x[0], apex3, apex4) && torus.contains(p); }, 4000000, rng);
    check(torus.body().solidVolume*regionShare23, torusInRegion, 5.0*torusRegionError/torusInRegion, "torus part in the throat region");

    // wetted surface: the cap of sphere 0 inside sphere 1, in the region
    std::normal_distribution<double> g(0.0, 1.0);
    const long surfaceSamples = 4000000;
    long capHits = 0;
    const double capAngle = std::acos(lens.radicalDistance/r[0]);
    const double capArea = 2.0*pi*r[0]*lens.capHeights[0];
    std::uniform_real_distribution<double> u(0.0, 1.0);
    for (long s = 0; s < surfaceSamples; ++s)
    {
        // uniform on the cap around the axis +x
        const double cosTheta = 1.0 - u(rng)*(1.0 - std::cos(capAngle));
        const double phi = 2.0*pi*u(rng);
        const double sinTheta = std::sqrt(1.0 - cosTheta*cosTheta);
        const P p{r[0]*cosTheta, r[0]*sinTheta*std::cos(phi), r[0]*sinTheta*std::sin(phi)};
        capHits += bipyramid(p, x[0], x[1], x[2], apex1, apex2);
    }
    const double capFraction = static_cast<double>(capHits)/surfaceSamples;
    const double capError = std::sqrt(capFraction*(1.0 - capFraction)/surfaceSamples);
    check(capFraction, regionShare01, 5.0*capError/regionShare01, "cap of sphere 0 in the throat region");
    check(-lensBody.wettedSurfaceChange, 2.0*pi*(r[0]*lens.capHeights[0] + r[1]*lens.capHeights[1]), 1e-14, "caps of the lens");
    check(capArea, 2.0*pi*r[0]*(r[0] - lens.radicalDistance), 1e-14, "cap area");
}

std::vector<P> latticeCenters(int n, double a)
{
    std::vector<P> centers;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
                centers.push_back(P{i*a, j*a, k*a});
    return centers;
}

//! Simple cubic lattice of spheres of radius R > a/2: closed forms with the lenses along the lattice edges
void testCoatedLattice()
{
    const int n = 5;
    const double a = 1.0, R = 0.52;
    const auto centers = latticeCenters(n, a);
    const std::vector<double> radii(centers.size(), R);
    const auto t = regularTriangulation(centers, radii);
    ExtractionRegion<double> region;
    region.lower = -0.1;
    region.upper = (n - 1)*a + 0.1;
    region.mergeDistance = 1e-9*a;
    const auto bodies = overlapCorrections(t, centers, radii);
    checkTrue(bodies.size() == std::size_t(3*n*n*(n - 1)), "one lens per lattice edge");
    const auto network = extractNetwork(t, centers, radii, region, bodies);

    const auto lens = sphereLens(R, R, a);
    const double h = R - 0.5*a;
    checkTrue(network.pores.size() == std::size_t((n - 1)*(n - 1)*(n - 1)), "one pore per unit cell");
    checkTrue(network.throats.size() == std::size_t(3*(n - 2)*(n - 1)*(n - 1)), "one throat per inner face");
    for (const auto& pore : network.pores)
        check(pore.volume, a*a*a - 4.0/3.0*pi*R*R*R + 3.0*lens.volume, 1e-12, "coated pore volume");

    const double fluidArea = a*a - pi*R*R + 2.0*lens.area;
    const double hydraulicRadius = (a*a*a/3.0 - 4.0*pi*R*R*R/9.0 + lens.volume)/(4.0*pi*R*R/3.0 - 4.0*pi*R*h);
    for (const auto& throat : network.throats)
    {
        check(throat.fluidArea, fluidArea, 1e-12, "coated throat fluid area");
        check(throat.hydraulicRadius, hydraulicRadius, 1e-12, "coated throat hydraulic radius");
    }
}

//! Random packing: without bodies and solid cells the network equals the one of extractNetwork
void testUnchanged(std::mt19937_64& rng)
{
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::vector<P> centers;
    std::vector<double> radii;
    const double size = 10.0;
    while (centers.size() < 300)
    {
        const double r = 0.4 + 0.2*u(rng);
        const P c{u(rng)*size, u(rng)*size, u(rng)*size};
        if (std::all_of(centers.begin(), centers.end(), [&, i = std::size_t(0)](const P& other) mutable {
                return (other - c).two_norm() >= r + radii[i++]; }))
        {
            centers.push_back(c);
            radii.push_back(r);
        }
    }
    const auto t = regularTriangulation(centers, radii);
    ExtractionRegion<double> region;
    region.lower = 1.5;
    region.upper = size - 1.5;
    const auto reference = extractNetwork(t, centers, radii, region);
    const auto network = extractNetwork(t, centers, radii, region, EdgeBodies<double>{});
    checkTrue(network.pores.size() == reference.pores.size() && network.throats.size() == reference.throats.size()
              && network.numClosedThroats == reference.numClosedThroats, "same numbers of pores and throats");
    bool equal = network.cellPores == reference.cellPores;
    for (std::size_t p = 0; equal && p < network.pores.size(); ++p)
        equal = network.pores[p].volume == reference.pores[p].volume && network.pores[p].label == reference.pores[p].label
                && network.pores[p].inscribedRadius == reference.pores[p].inscribedRadius;
    for (std::size_t k = 0; equal && k < network.throats.size(); ++k)
    {
        const auto& t0 = network.throats[k];
        const auto& t1 = reference.throats[k];
        equal = t0.pores == t1.pores && t0.fluidArea == t1.fluidArea && t0.hydraulicRadius == t1.hydraulicRadius
                && t0.length == t1.length && t0.shapeFactor == t1.shapeFactor;
    }
    checkTrue(equal, "pores and throats identical to extractNetwork");
}

//! Solid unit cells of a lattice: a solid layer splits the network, an enclosed pair of cells is pruned
void testSolidCells()
{
    const int n = 7;
    const double a = 1.0, r = 0.45;
    const auto centers = latticeCenters(n, a);
    const std::vector<double> radii(centers.size(), r);
    const auto t = regularTriangulation(centers, radii);
    // the outer unit cells reach beyond the region and are boundary pores
    ExtractionRegion<double> region;
    region.lower = 0.2*a;
    region.upper = (n - 1)*a - 0.2*a;
    region.mergeDistance = 1e-9*a;

    using Cube = std::array<int, 3>;
    std::vector<Cube> solid{{1, 1, 2}, {1, 2, 1}, {1, 2, 3}, {0, 2, 2}, {2, 2, 2},
                            {1, 4, 2}, {1, 3, 1}, {1, 3, 3}, {0, 3, 2}, {2, 3, 2}};
    for (int j = 0; j < n - 1; ++j)
        for (int k = 0; k < n - 1; ++k)
            solid.push_back({3, j, k});
    std::vector<char> solidCells(t.cells.size(), 0);
    for (std::size_t c = 0; c < t.cells.size(); ++c)
    {
        P centroid(0.0);
        for (const int v : t.cells[c])
            centroid += centers[v];
        centroid /= 4.0;
        const Cube cube{static_cast<int>(std::floor(centroid[0])), static_cast<int>(std::floor(centroid[1])),
                        static_cast<int>(std::floor(centroid[2]))};
        solidCells[c] = std::find(solid.begin(), solid.end(), cube) != solid.end();
    }

    const auto network = extractNetwork(t, centers, radii, region, EdgeBodies<double>{}, solidCells);
    const std::size_t numCubes = (n - 1)*(n - 1)*(n - 1);
    checkTrue(network.pores.size() == numCubes - solid.size(), "solid cells are no pores");
    const auto connected = connectedNetwork(network, {region.boundaryLabels[0], region.boundaryLabels[1]});
    checkTrue(connected.pores.size() == network.pores.size() - 2, "enclosed pair of pores removed");
    checkTrue(connected.throats.size() == network.throats.size() - 1, "throat of the enclosed pair removed");
    bool valid = true;
    for (const auto& throat : connected.throats)
        valid = valid && throat.pores[0] < connected.pores.size() && throat.pores[1] < connected.pores.size();
    for (const auto p : connected.cellPores)
        valid = valid && p < static_cast<long>(connected.pores.size());
    checkTrue(valid, "renumbered throats and cells");
}

void testTripleOverlap()
{
    const std::vector<P> centers{P{0.0, 0.0, 0.0}, P{1.0, 0.0, 0.0}, P{0.5, std::sqrt(3.0)/2.0, 0.0},
                                 P{0.5, std::sqrt(3.0)/6.0, std::sqrt(2.0/3.0)}};
    const auto t = regularTriangulation(centers, std::vector<double>(4, 0.5));
    bool thrown = false;
    try { overlapCorrections(t, centers, std::vector<double>(4, 0.6)); }
    catch (const Dune::InvalidStateException&) { thrown = true; }
    checkTrue(thrown, "three overlapping spheres are rejected");
    overlapCorrections(t, centers, std::vector<double>(4, 0.55));

    const std::array<P, 3> x{centers[0], centers[1], centers[2]};
    const std::array<double, 3> radius{0.3, 0.5, 0.4};
    const auto f = facetPowerCenter(x, radius);
    const double p0 = (f - x[0]).two_norm2() - radius[0]*radius[0];
    check((f - x[1]).two_norm2() - radius[1]*radius[1], p0, 1e-13, "facet power centre, sphere 1");
    check((f - x[2]).two_norm2() - radius[2]*radius[2], p0, 1e-13, "facet power centre, sphere 2");
    checkTrue(std::abs(f[2]) < 1e-15, "facet power centre in the plane");
}

} // end namespace

int main(int argc, char** argv)
{
    Dumux::initialize(argc, argv);
    std::mt19937_64 rng(42);
    testLens(rng);
    testShares(rng);
    testCoatedLattice();
    testUnchanged(rng);
    testSolidCells();
    testTripleOverlap();
    if (numFailures > 0)
    {
        std::cout << numFailures << " checks failed\n";
        return 1;
    }
    std::cout << "all checks passed\n";
    return 0;
}
