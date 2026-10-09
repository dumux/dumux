// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkModels
 * \brief Test of the pore network extracted from the regular triangulation of a sphere packing
 */
#include <config.h>

#include <algorithm>
#include <array>
#include <functional>
#include <cmath>
#include <iostream>
#include <map>
#include <numeric>
#include <numbers>
#include <random>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/io/grid/porenetwork/gridmanager.hh>
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

struct Packing
{
    std::vector<P> centers;
    std::vector<double> radii;
    P lower, upper;
};

// random sequential addition of non-overlapping spheres, largest first
Packing randomPacking(int n, double solidFraction, unsigned int seed)
{
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> radius(0.4, 0.6), unit(0.0, 1.0);
    Packing p;
    for (int i = 0; i < n; ++i)
        p.radii.push_back(radius(gen));
    std::sort(p.radii.begin(), p.radii.end(), std::greater<double>{});
    double solid = 0.0;
    for (const double r : p.radii)
        solid += 4.0/3.0*pi*r*r*r;
    const double length = std::cbrt(solid/solidFraction);
    p.lower = 0.0;
    p.upper = length;
    for (int i = 0; i < n; ++i)
    {
        const double r = p.radii[i];
        for (int attempt = 0; ; ++attempt)
        {
            if (attempt > 100000)
                DUNE_THROW(Dune::Exception, "random packing failed");
            const P x{r + (length - 2*r)*unit(gen), r + (length - 2*r)*unit(gen), r + (length - 2*r)*unit(gen)};
            bool free = true;
            for (std::size_t j = 0; j < p.centers.size() && free; ++j)
                free = (x - p.centers[j]).two_norm() >= r + p.radii[j];
            if (free)
            {
                p.centers.push_back(x);
                break;
            }
        }
    }
    return p;
}

std::array<P, 4> cellPoints(const Triangulation& t, const Packing& p, std::size_t c)
{
    std::array<P, 4> x;
    for (int j = 0; j < 4; ++j)
        x[j] = p.centers[t.cells[c][j]];
    return x;
}

void testTriangulation()
{
    const auto packing = randomPacking(500, 0.25, 1);
    const auto t = regularTriangulation(packing.centers, packing.radii);
    checkTrue(t.numHiddenSpheres == 0, "no hidden spheres in a packing without overlaps");

    // solid angles around every sphere not on the convex hull sum to 4 pi
    std::vector<double> angleSum(packing.centers.size(), 0.0);
    std::vector<bool> onHull(packing.centers.size(), false);
    double cellVolume = 0.0, hullVolume = 0.0;
    for (std::size_t c = 0; c < t.cells.size(); ++c)
    {
        const auto x = cellPoints(t, packing, c);
        cellVolume += tetrahedronVolume(x);
        for (int j = 0; j < 4; ++j)
        {
            angleSum[t.cells[c][j]] += solidAngle(x[j], x[(j+1)%4], x[(j+2)%4], x[(j+3)%4]);
            if (t.neighbors[c][j] < 0)
            {
                // hull facet opposite to vertex j, oriented away from it
                const auto& a = x[(j+1)%4];
                const auto& b = x[(j+2)%4];
                const auto& d = x[(j+3)%4];
                auto normal = Dumux::PoreNetwork::SpherePacking::Detail::cross(b - a, d - a);
                if (normal*(x[j] - a) > 0.0)
                    normal *= -1.0;
                hullVolume += (a*normal)/6.0;
                for (int k = 1; k < 4; ++k)
                    onHull[t.cells[c][(j+k)%4]] = true;
            }
        }
    }
    int numInterior = 0;
    for (std::size_t i = 0; i < angleSum.size(); ++i)
    {
        if (onHull[i])
            continue;
        ++numInterior;
        check(angleSum[i], 4.0*pi, 1e-12, "solid angles around an interior sphere");
    }
    checkTrue(numInterior > 100, "enough interior spheres");
    check(cellVolume, hullVolume, 1e-12, "tetrahedra tile the convex hull");
}

void testCubicLattice()
{
    // simple cubic lattice: all tetrahedra of a unit cell share the power centre and merge into one pore
    const int n = 5;
    const double a = 1.0, r = 0.45;
    Packing p;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
            {
                p.centers.push_back(P{i*a, j*a, k*a});
                p.radii.push_back(r);
            }

    const auto t = regularTriangulation(p.centers, p.radii);
    ExtractionRegion<double> region;
    region.lower = -0.1;
    region.upper = (n - 1)*a + 0.1;
    region.mergeDistance = 1e-9*a;
    const auto network = extractNetwork(t, p.centers, p.radii, region);

    const int numCells = (n - 1)*(n - 1)*(n - 1);
    checkTrue(network.pores.size() == std::size_t(numCells), "one pore per unit cell");
    checkTrue(network.throats.size() == std::size_t(3*(n - 2)*(n - 1)*(n - 1)), "one throat per inner face");
    checkTrue(network.numClosedThroats == 0 && network.numIsolatedPores == 0, "no closed throats and isolated pores");

    for (const auto& pore : network.pores)
    {
        check(pore.volume, a*a*a - 4.0/3.0*pi*r*r*r, 1e-12, "pore volume of the unit cell");
        check(pore.bulkVolume, a*a*a, 1e-12, "bulk volume of the unit cell");
        check(pore.inscribedRadius, std::sqrt(3.0)/2.0*a - r, 1e-12, "inscribed sphere of the unit cell");
        checkTrue(pore.label == -1, "interior label");
        for (int axis = 0; axis < 3; ++axis)
        {
            const double offset = pore.position[axis]/a - std::floor(pore.position[axis]/a);
            check(offset, 0.5, 1e-12, "pore at the centre of the unit cell");
        }
    }

    const double hydraulicRadius = (a*a*a/3.0 - 4.0*pi*r*r*r/9.0)/(4.0*pi*r*r/3.0);
    for (const auto& throat : network.throats)
    {
        check(throat.length, a, 1e-12, "throat length");
        check(throat.fluidArea, a*a - pi*r*r, 1e-12, "fluid area of a face");
        check(throat.inscribedRadius, std::sqrt(2.0)/2.0*a - r, 1e-12, "constriction of a face");
        check(throat.hydraulicRadius, hydraulicRadius, 1e-12, "hydraulic radius of a face");
    }

    // linear pressure: every face weighs a^2/4 per corner, so an interior sphere carries the pressure
    // gradient times the volume a^3 per sphere
    const P gradient{0.3, -1.1, 0.7};
    std::vector<double> pressure;
    for (const auto& pore : network.pores)
        pressure.push_back(gradient*pore.position);
    const auto forces = fluidForces(network, pressure, p.centers.size());
    for (std::size_t s = 0; s < p.centers.size(); ++s)
    {
        const auto& x = p.centers[s];
        bool interior = true;
        for (int axis = 0; axis < 3; ++axis)
            interior = interior && x[axis] > 0.5*a && x[axis] < (n - 1.5)*a;
        if (interior)
            for (int axis = 0; axis < 3; ++axis)
                check(forces[s][axis], -a*a*a*gradient[axis], 1e-12, "force of a linear pressure on an interior sphere");
    }
}

void testRegionAndDgf()
{
    const auto packing = randomPacking(2000, 0.25, 2);
    const auto t = regularTriangulation(packing.centers, packing.radii);
    const double margin = 2.0;
    ExtractionRegion<double> region;
    region.lower = packing.lower + P(margin);
    region.upper = packing.upper - P(margin);
    const auto network = extractNetwork(t, packing.centers, packing.radii, region);
    std::cout << "random packing: " << network.pores.size() << " pores, " << network.throats.size() << " throats, "
              << network.numClosedThroats << " closed throats, " << network.numIsolatedPores << " isolated pores\n";

    std::map<int, int> labelCount;
    for (const auto& pore : network.pores)
    {
        checkTrue(pore.inscribedRadius > 0.0 && pore.volume > 0.0, "positive pore radius and volume");
        ++labelCount[pore.label];
        if (pore.label > 0)
        {
            const int face = pore.label - 1;
            const int axis = face/2;
            const double distance = face % 2 == 0 ? pore.position[axis] - region.lower[axis] : region.upper[axis] - pore.position[axis];
            checkTrue(std::abs(distance) < 2.0, "boundary pore close to its face");
        }
    }
    for (int label = 1; label <= 6; ++label)
        checkTrue(labelCount[label] > 0, "pores on every boundary face");
    for (const auto& throat : network.throats)
        checkTrue(throat.length > 0.0 && throat.fluidArea > 0.0 && throat.inscribedRadius > 0.0
                  && throat.hydraulicRadius > 0.0 && throat.shapeFactor > 0.0, "positive throat parameters");

    // pressure forces: none for a uniform pressure, unit normals of the facets pointing into the first
    // pore, weights summing to the facet area if no circle crosses an edge
    const auto uniformForces = fluidForces(network, std::vector<double>(network.pores.size(), 3.7), packing.centers.size());
    checkTrue(std::all_of(uniformForces.begin(), uniformForces.end(), [](const P& f) { return f.two_norm() == 0.0; }),
              "no force for a uniform pressure");
    checkTrue(network.facets.size() >= network.throats.size(), "a facet for every throat");
    for (const auto& f : network.facets)
    {
        std::array<P, 3> x;
        std::array<double, 3> r;
        for (int k = 0; k < 3; ++k)
        {
            x[k] = packing.centers[f.vertices[k]];
            r[k] = packing.radii[f.vertices[k]];
        }
        check(f.normal.two_norm(), 1.0, 1e-14, "unit facet normal");
        check(1.0 + f.normal*(x[1] - x[0])/(x[1] - x[0]).two_norm(), 1.0, 1e-13, "facet normal orthogonal to the facet");
        checkTrue(f.normal*(network.pores[f.pores[0]].position - network.pores[f.pores[1]].position) > 0.0,
                  "facet normal pointing into the first pore");
        bool crossing = false;
        for (int k = 0; k < 3; ++k)
            crossing = crossing || circularSegmentBeyondEdge(x[k], r[k], x[(k+1)%3], x[(k+2)%3]) > 0.0;
        if (!crossing)
            check(f.forceWeights[0] + f.forceWeights[1] + f.forceWeights[2], triangleArea(x[0], x[1], x[2]),
                  1e-12, "force weights sum to the facet area");
    }

    // round trip through the DGF reader of the pore-network grid manager
    const std::string fileName = "spherepacking_network.dgf";
    writeDgf(fileName, network);
    Dumux::Parameters::init([&](auto& tree) {
        tree["Grid.File"] = fileName;
        tree["Grid.Sanitize"] = "false";
    });
    Dumux::PoreNetwork::GridManager<3> gridManager;
    gridManager.init();
    const auto gridView = gridManager.grid().leafGridView();
    const auto gridData = gridManager.getGridData();
    checkTrue(gridView.size(1) == network.pores.size(), "number of pores read back");
    checkTrue(gridView.size(0) == network.throats.size(), "number of throats read back");

    std::map<std::array<long, 3>, std::size_t> poreAt;
    const auto key = [](const P& x) {
        return std::array<long, 3>{std::lround(x[0]*1e9), std::lround(x[1]*1e9), std::lround(x[2]*1e9)};
    };
    for (std::size_t i = 0; i < network.pores.size(); ++i)
        poreAt[key(network.pores[i].position)] = i;
    for (const auto& vertex : vertices(gridView))
    {
        const auto& pore = network.pores.at(poreAt.at(key(vertex.geometry().center())));
        check(gridData->getParameter(vertex, "PoreVolume"), pore.volume, 1e-15, "pore volume read back");
        check(gridData->getParameter(vertex, "PoreInscribedRadius"), pore.inscribedRadius, 1e-15, "pore radius read back");
        checkTrue(gridData->getParameter(vertex, "PoreLabel") == pore.label, "pore label read back");
    }
    double length = 0.0, area = 0.0, hydraulic = 0.0, networkLength = 0.0, networkArea = 0.0, networkHydraulic = 0.0;
    for (const auto& element : elements(gridView))
    {
        length += gridData->getParameter(element, "ThroatLength");
        area += gridData->getParameter(element, "ThroatCrossSectionalArea");
        hydraulic += gridData->getParameter(element, "ThroatHydraulicRadius");
    }
    for (const auto& throat : network.throats)
    {
        networkLength += throat.length;
        networkArea += throat.fluidArea;
        networkHydraulic += throat.hydraulicRadius;
    }
    check(length, networkLength, 1e-13, "throat lengths read back");
    check(area, networkArea, 1e-13, "throat areas read back");
    check(hydraulic, networkHydraulic, 1e-13, "hydraulic radii read back");
}

void testWalls()
{
    const auto packing = randomPacking(1000, 0.25, 4);
    Walls<double> walls;
    walls.lower = packing.lower;
    walls.upper = packing.upper;
    walls.numSpheres = packing.centers.size();
    const auto t = regularTriangulation(packing.centers, packing.radii, walls);
    checkTrue(t.numHiddenSpheres == 0, "no hidden spheres with walls");
    const auto network = extractNetwork(t, packing.centers, packing.radii, walls);
    std::cout << "walls: " << network.pores.size() << " pores, " << network.throats.size() << " throats, "
              << network.numClosedThroats << " closed throats, " << network.numIsolatedPores << " isolated pores, "
              << network.wallPressures.size() << " wall pressure terms\n";

    const auto extent = packing.upper - packing.lower;
    const double boxVolume = extent[0]*extent[1]*extent[2];
    double bulk = 0.0, solid = 0.0, sphereVolume = 0.0;
    for (const auto& pore : network.pores)
    {
        bulk += pore.bulkVolume;
        solid += pore.bulkVolume - pore.volume;
    }
    for (const double r : packing.radii)
        sphereVolume += 4.0/3.0*pi*r*r*r;
    check(bulk, boxVolume, 1e-10, "tetrahedra with walls tile the box");
    check(solid, sphereVolume, 1e-10, "solid parts of the pores add up to the spheres");

    std::array<int, 6> wallPores{};
    for (const auto& pore : network.pores)
        for (int k = 0; k < 6; ++k)
            wallPores[k] += pore.label > 0 && (pore.label & (1 << k));
    for (int k = 0; k < 6; ++k)
        checkTrue(wallPores[k] > 0, "pores on every wall");
    for (const auto& pore : network.pores)
        checkTrue(pore.inscribedRadius > 0.0 && pore.volume > 0.0, "positive pore radius and volume with walls");
    for (const auto& throat : network.throats)
        checkTrue(throat.length > 0.0 && throat.fluidArea > 0.0 && throat.inscribedRadius > 0.0
                  && throat.hydraulicRadius > 0.0 && throat.shapeFactor > 0.0, "positive throat parameters with walls");

    // after displacing the spheres and moving the walls, the tetrahedra of the triangulation still tile the box
    std::mt19937 gen(9);
    std::uniform_real_distribution<double> shift(-1e-4, 1e-4);
    auto moved = packing.centers;
    for (auto& x : moved)
        for (int c = 0; c < 3; ++c)
            x[c] += shift(gen);
    auto movedWalls = walls;
    movedWalls.lower += P{-1e-3*extent[0], 2e-3*extent[1], 0.5e-3*extent[2]};
    movedWalls.upper += P{2e-3*extent[0], -1e-3*extent[1], 1e-3*extent[2]};
    const auto movedExtent = movedWalls.upper - movedWalls.lower;
    const auto movedVolumes = poreBulkVolumes(network, t, moved, packing.radii, &movedWalls);
    check(std::accumulate(movedVolumes.begin(), movedVolumes.end(), 0.0), movedExtent[0]*movedExtent[1]*movedExtent[2],
          1e-10, "bulk volumes after a deformation tile the moved box");
    const auto unmovedVolumes = poreBulkVolumes(network, t, packing.centers, packing.radii, &walls);
    for (std::size_t i = 0; i < network.pores.size(); ++i)
        check(unmovedVolumes[i], network.pores[i].bulkVolume, 1e-14, "bulk volumes recomputed without deformation");

    // a uniform pressure leaves the spheres free of force and pushes every wall outwards with p times its area
    const double p = 2.5;
    const auto forces = fluidForces(network, std::vector<double>(network.pores.size(), p), packing.centers.size() + 6);
    double sphereForce = 0.0;
    for (std::size_t s = 0; s < packing.centers.size(); ++s)
        sphereForce = std::max(sphereForce, forces[s].two_norm());
    checkTrue(sphereForce == 0.0, "no force on the spheres for a uniform pressure with walls");
    for (int k = 0; k < 6; ++k)
    {
        const int axis = k/2;
        const double area = extent[(axis + 1)%3]*extent[(axis + 2)%3];
        const auto error = forces[packing.centers.size() + k] + p*area*walls.inwardNormal(k);
        check(1.0 + error.two_norm()/(p*area), 1.0, 1e-12, "wall force of a uniform pressure");
    }
}

} // end anonymous namespace

int main(int argc, char** argv)
{
    Dumux::initialize(argc, argv);

    testTriangulation();
    testCubicLattice();
    testRegionAndDgf();
    testWalls();

    if (numFailures > 0)
        DUNE_THROW(Dune::Exception, numFailures << " checks failed");
    std::cout << "all checks passed\n";
    return 0;
}
