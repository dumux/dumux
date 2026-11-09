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
#include <queue>
#include <random>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/io/grid/porenetwork/dgfwriter.hh>
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

    // the throat regions of the interior facets and the pyramids from the power centres over the hull
    // facets, negative for a power centre beyond its facet, partition the void volume and the wetted
    // surface of the tetrahedra
    double voidVolume = 0.0, wettedSurface = 0.0, regionVoid = 0.0, regionSurface = 0.0, scale = 0.0;
    std::vector<P> powerCenters;
    for (std::size_t c = 0; c < t.cells.size(); ++c)
    {
        const auto x = cellPoints(t, packing, c);
        std::array<double, 4> r;
        for (int j = 0; j < 4; ++j)
            r[j] = packing.radii[t.cells[c][j]];
        powerCenters.push_back(powerCenter(x, r));
        voidVolume += tetrahedronVoidVolume(x, r);
        for (int j = 0; j < 4; ++j)
            wettedSurface += solidAngle(x[j], x[(j+1)%4], x[(j+2)%4], x[(j+3)%4])*r[j]*r[j];
    }
    int numFootOutside = 0;
    for (std::size_t c = 0; c < t.cells.size(); ++c)
        for (int j = 0; j < 4; ++j)
        {
            const int n = t.neighbors[c][j];
            if (n >= 0 && static_cast<std::size_t>(n) < c)
                continue;
            std::array<P, 3> x;
            std::array<double, 3> r;
            for (int k = 0, m = 0; k < 4; ++k)
                if (k != j)
                {
                    x[m] = packing.centers[t.cells[c][k]];
                    r[m++] = packing.radii[t.cells[c][k]];
                }
            const auto normal = Dumux::PoreNetwork::SpherePacking::Detail::cross(x[1] - x[0], x[2] - x[0]);
            if (n >= 0)
            {
                const auto region = throatRegion(x, r, powerCenters[c], powerCenters[n]);
                regionVoid += region.voidVolume;
                regionSurface += region.solidSurface;
                scale += std::abs(region.volume);
                const auto foot = powerCenters[c] - ((powerCenters[c] - x[0])*normal/normal.two_norm2())*normal;
                bool inside = true;
                for (int k = 0; k < 3; ++k)
                    inside = inside && Dumux::PoreNetwork::SpherePacking::Detail::cross(x[(k+1)%3] - x[k], foot - x[k])*normal >= 0.0;
                numFootOutside += !inside;
                continue;
            }
            const auto& opposite = packing.centers[t.cells[c][j]];
            const double sign = (normal*(opposite - x[0]))*(normal*(powerCenters[c] - x[0])) >= 0.0 ? 1.0 : -1.0;
            const double volume = tetrahedronVolume(std::array<P, 4>{x[0], x[1], x[2], powerCenters[c]});
            scale += volume;
            regionVoid += sign*volume;
            for (int k = 0; k < 3; ++k)
            {
                const double omega = solidAngle(x[k], x[(k+1)%3], x[(k+2)%3], powerCenters[c]);
                regionVoid -= sign*omega*r[k]*r[k]*r[k]/3.0;
                regionSurface += sign*omega*r[k]*r[k];
            }
        }
    std::cout << "throat regions: " << numFootOutside << " interior facets whose dual edge misses the facet, relative "
              << "deviation of the void volume " << regionVoid/voidVolume - 1.0 << " and of the wetted surface "
              << regionSurface/wettedSurface - 1.0 << ", summed region volumes " << scale/voidVolume << " times the void volume\n";
    checkTrue(numFootOutside > 100, "enough facets whose dual edge misses the facet");
    check(regionVoid, voidVolume, 1e-12*scale/voidVolume, "throat regions partition the void volume");
    check(regionSurface, wettedSurface, 1e-12*scale/voidVolume, "throat regions partition the wetted surface");
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

    // every facet joins the tetrahedra of its two pores, which share its three vertices
    checkTrue(network.facets.size() >= network.throats.size(), "a facet for every throat");
    for (const auto& f : network.facets)
        for (int side = 0; side < 2; ++side)
        {
            checkTrue(network.cellPores[f.cells[side]] == static_cast<long>(f.pores[side]), "facet tetrahedron in its pore");
            const auto& cell = t.cells[f.cells[side]];
            for (const int v : f.vertices)
                checkTrue(std::find(cell.begin(), cell.end(), v) != cell.end(), "facet vertex in its tetrahedron");
        }

    // grid from the network data, written by the DGF writer of the pore-network grids and read back
    Dumux::Parameters::init([](auto& tree) { tree["Grid.Sanitize"] = "false"; });
    const std::string fileName = "spherepacking_network.dgf";
    {
        Dumux::PoreNetwork::GridManager<3> dataGridManager;
        dataGridManager.init(networkData(network));
        Dumux::PoreNetwork::writeDgf(fileName, dataGridManager.grid().leafGridView(), *dataGridManager.getGridData());
    }
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

} // end anonymous namespace

int main(int argc, char** argv)
{
    Dumux::initialize(argc, argv);

    testTriangulation();
    testCubicLattice();
    testRegionAndDgf();

    if (numFailures > 0)
        DUNE_THROW(Dune::Exception, numFailures << " checks failed");
    std::cout << "all checks passed\n";
    return 0;
}
