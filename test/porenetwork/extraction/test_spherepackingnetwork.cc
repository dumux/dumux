// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkModels
 * \brief Test of the pore network of a sphere packing and of single-phase flow through it
 *
 * Simple cubic lattice: the pores and throats of the unit cells and the permeability A R_h^2/(2 a^2) with the
 * fluid area A and hydraulic radius R_h of a face and the lattice constant a. Random packing: the tetrahedra tile
 * the convex hull, the throat regions partition the void volume and the wetted surface of the tetrahedra, inflow
 * equals outflow and the flux of the model equals the flux from the throat transmissibilities, for a grid from the
 * network data as for the same grid read from a DGF file, and for given transmissibilities.
 */
#include <config.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <numbers>
#include <random>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/io/grid/porenetwork/dgfwriter.hh>
#include <dumux/io/grid/porenetwork/gridmanager.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/porenetwork/common/boundaryflux.hh>
#include <dumux/porenetwork/common/labels.hh>
#include <dumux/porenetwork/extraction/spherepackingtriangulation.hh>

#include "properties.hh"

namespace {

using namespace Dumux::PoreNetwork::SpherePacking;
using P = Point<double>;
using Labels = Dumux::PoreNetwork::Labels;
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
    double length;
};

// random sequential addition of non-overlapping spheres in a cube, largest first
Packing randomPacking(int n, double solidFraction, unsigned int seed)
{
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> radius(0.4e-3, 0.6e-3), unit(0.0, 1.0);
    Packing p;
    for (int i = 0; i < n; ++i)
        p.radii.push_back(radius(gen));
    std::sort(p.radii.begin(), p.radii.end(), std::greater<double>{});
    double solid = 0.0;
    for (const double r : p.radii)
        solid += 4.0/3.0*pi*r*r*r;
    p.length = std::cbrt(solid/solidFraction);
    for (const double r : p.radii)
        while (true)
        {
            const P x{r + (p.length - 2*r)*unit(gen), r + (p.length - 2*r)*unit(gen), r + (p.length - 2*r)*unit(gen)};
            bool free = true;
            for (std::size_t j = 0; j < p.centers.size() && free; ++j)
                free = (x - p.centers[j]).two_norm() >= r + p.radii[j];
            if (free)
            {
                p.centers.push_back(x);
                break;
            }
        }
    return p;
}

struct Flow
{
    double inflow, outflow; //!< volume fluxes through the inlet and outlet pores computed by the model
    double networkFlux; //!< volume flux out of the inlet pores from the throat transmissibilities and the pressures
};

//! flow from the inlet pores (pressure 1) to the outlet pores (pressure 0) with the single-phase pore-network model
template<class TypeTag = Dumux::Properties::TTag::SpherePackingPermeability>
Flow solveFlow(Dumux::PoreNetwork::GridManager<3>& gridManager, const Network<double>& network)
{
    using namespace Dumux;
    const auto& gridView = gridManager.grid().leafGridView();
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(gridView, *gridManager.getGridData());
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;
    auto spatialParams = std::make_shared<SpatialParams>(gridGeometry, *gridManager.getGridData());
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry, spatialParams);

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    x = 0.0;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);
    using Assembler = FVAssembler<TypeTag, DiffMethod::analytic>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    auto A = std::make_shared<JacobianMatrix>();
    auto r = std::make_shared<SolutionVector>();
    assembler->setLinearSystem(A, r);
    assembler->assembleJacobianAndResidual(x);
    (*r) *= -1.0;
    UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>().solve(*A, x, *r);
    gridVariables->update(x);

    // the boundary flux is positive for flow out of the network
    const double density = getParam<double>("Component.LiquidDensity");
    const auto boundaryFlux = PoreNetwork::BoundaryFlux(*gridVariables, assembler->localResidual(), x);
    Flow flow{-boundaryFlux.getFlux(std::vector<int>{Labels::inlet}).totalFlux[0]/density,
              boundaryFlux.getFlux(std::vector<int>{Labels::outlet}).totalFlux[0]/density, 0.0};

    std::map<std::array<long, 3>, double> pressureAt;
    const auto key = [](const P& p) { return std::array<long, 3>{std::lround(p[0]*1e12), std::lround(p[1]*1e12), std::lround(p[2]*1e12)}; };
    for (const auto& vertex : vertices(gridView))
        pressureAt[key(vertex.geometry().center())] = x[gridGeometry->vertexMapper().index(vertex)][0];
    for (const auto& throat : network.throats)
        for (int side = 0; side < 2; ++side)
        {
            const auto& pore = network.pores[throat.pores[side]];
            const auto& other = network.pores[throat.pores[1-side]];
            if (pore.label == Labels::inlet && other.label != Labels::inlet)
                flow.networkFlux += throatTransmissibility(throat)/getParam<double>("Component.LiquidDynamicViscosity")
                                    *(pressureAt.at(key(pore.position)) - pressureAt.at(key(other.position)));
        }
    return flow;
}

template<class TypeTag = Dumux::Properties::TTag::SpherePackingPermeability>
Flow solveFlow(const Network<double>& network)
{
    Dumux::PoreNetwork::GridManager<3> gridManager;
    gridManager.init(networkData(network));
    return solveFlow<TypeTag>(gridManager, network);
}

//! the x faces are inlet and outlet, the other faces impermeable
ExtractionRegion<double> flowRegion(const P& lower, const P& upper)
{ return {lower, upper, {Labels::inlet, Labels::outlet, Labels::noflow, Labels::noflow, Labels::noflow, Labels::noflow}}; }

void testCubicLattice()
{
    // all tetrahedra of a unit cell share the power centre and merge into one pore
    const int n = 6;
    const double a = 1e-3, r = 0.45e-3;
    std::vector<P> centers;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
                centers.push_back(P{i*a, j*a, k*a});
    const std::vector<double> radii(centers.size(), r);
    auto region = flowRegion(P(-0.1*a), P((n - 1 + 0.1)*a));
    region.mergeDistance = 1e-9*a;
    auto network = extractNetwork(regularTriangulation(centers, radii), centers, radii, region);

    checkTrue(network.pores.size() == std::size_t((n - 1)*(n - 1)*(n - 1)), "lattice: one pore per unit cell");
    checkTrue(network.throats.size() == std::size_t(3*(n - 2)*(n - 1)*(n - 1)), "lattice: one throat per inner face");
    for (const auto& pore : network.pores)
    {
        check(pore.volume, a*a*a - 4.0/3.0*pi*r*r*r, 1e-12, "lattice: pore volume of the unit cell");
        check(pore.inscribedRadius, std::sqrt(3.0)/2.0*a - r, 1e-12, "lattice: inscribed sphere of the unit cell");
        for (int axis = 0; axis < 3; ++axis)
            check(pore.position[axis]/a - std::floor(pore.position[axis]/a), 0.5, 1e-12, "lattice: pore at the centre of the unit cell");
    }
    const double fluidArea = a*a - pi*r*r;
    const double hydraulicRadius = (a*a*a/3.0 - 4.0*pi*r*r*r/9.0)/(4.0*pi*r*r/3.0);
    for (const auto& throat : network.throats)
    {
        check(throat.length, a, 1e-12, "lattice: throat length");
        check(throat.fluidArea, fluidArea, 1e-12, "lattice: fluid area of a face");
        check(throat.inscribedRadius, std::sqrt(2.0)/2.0*a - r, 1e-12, "lattice: constriction of a face");
        check(throat.hydraulicRadius, hydraulicRadius, 1e-12, "lattice: hydraulic radius of a face");
    }

    // inlet and outlet: the first and last layer of unit cells, a distance (n - 2) a apart
    for (auto& pore : network.pores)
        pore.label = pore.position[0] < a ? Labels::inlet : (pore.position[0] > (n - 2)*a ? Labels::outlet : Labels::interior);
    const auto flow = solveFlow(network);
    const double permeability = Dumux::getParam<double>("Component.LiquidDynamicViscosity")*flow.inflow*(n - 2)*a/((n - 1)*a*(n - 1)*a);
    check(permeability, fluidArea*hydraulicRadius*hydraulicRadius/(2.0*a*a), 1e-10, "lattice: permeability");
    check(flow.outflow, flow.inflow, 1e-10, "lattice: outflow equals inflow");
}

void testRandomPacking()
{
    const auto [centers, radii, length] = randomPacking(2000, 0.25, 3);
    const auto t = regularTriangulation(centers, radii);
    checkTrue(t.numHiddenSpheres == 0, "random packing: no hidden spheres without overlaps");

    // the tetrahedra tile the convex hull; the throat regions of the interior facets and the pyramids from the
    // power centres over the hull facets, negative for a power centre beyond its facet, partition the void volume
    // and the wetted surface of the tetrahedra
    double cellVolume = 0.0, hullVolume = 0.0, voidVolume = 0.0, wettedSurface = 0.0, regionVoid = 0.0, regionSurface = 0.0;
    std::vector<P> powerCenters;
    for (const auto& cell : t.cells)
    {
        const std::array<P, 4> x{centers[cell[0]], centers[cell[1]], centers[cell[2]], centers[cell[3]]};
        const std::array<double, 4> r{radii[cell[0]], radii[cell[1]], radii[cell[2]], radii[cell[3]]};
        powerCenters.push_back(powerCenter(x, r));
        cellVolume += tetrahedronVolume(x);
        voidVolume += tetrahedronVoidVolume(x, r);
        for (int j = 0; j < 4; ++j)
            wettedSurface += solidAngle(x[j], x[(j+1)%4], x[(j+2)%4], x[(j+3)%4])*r[j]*r[j];
    }
    for (std::size_t c = 0; c < t.cells.size(); ++c)
        for (int j = 0; j < 4; ++j)
        {
            const int n = t.neighbors[c][j];
            if (n >= 0 && static_cast<std::size_t>(n) < c)
                continue;
            std::array<P, 3> x;
            std::array<double, 3> r;
            for (int k = 1; k < 4; ++k)
            {
                x[k-1] = centers[t.cells[c][(j+k)%4]];
                r[k-1] = radii[t.cells[c][(j+k)%4]];
            }
            if (n >= 0)
            {
                const auto region = throatRegion(x, r, powerCenters[c], powerCenters[n]);
                regionVoid += region.voidVolume;
                regionSurface += region.solidSurface;
                continue;
            }
            auto outward = Dumux::crossProduct(x[1] - x[0], x[2] - x[0]);
            if (outward*(centers[t.cells[c][j]] - x[0]) > 0.0)
                outward *= -1.0;
            hullVolume += (x[0]*outward)/6.0;
            const double sign = outward*(powerCenters[c] - x[0]) <= 0.0 ? 1.0 : -1.0;
            regionVoid += sign*tetrahedronVolume(std::array<P, 4>{x[0], x[1], x[2], powerCenters[c]});
            for (int k = 0; k < 3; ++k)
            {
                const double omega = solidAngle(x[k], x[(k+1)%3], x[(k+2)%3], powerCenters[c]);
                regionVoid -= sign*omega*r[k]*r[k]*r[k]/3.0;
                regionSurface += sign*omega*r[k]*r[k];
            }
        }
    check(cellVolume, hullVolume, 1e-12, "random packing: tetrahedra tile the convex hull");
    check(regionVoid, voidVolume, 1e-10, "random packing: throat regions partition the void volume");
    check(regionSurface, wettedSurface, 1e-10, "random packing: throat regions partition the wetted surface");

    const auto network = extractNetwork(t, centers, radii, flowRegion(P(2e-3), P(length - 2e-3)));
    for (const auto& throat : network.throats)
        checkTrue(throat.length > 0.0 && throat.fluidArea > 0.0 && throat.inscribedRadius > 0.0
                  && throat.hydraulicRadius > 0.0 && throat.shapeFactor > 0.0, "random packing: positive throat parameters");
    for (const auto& f : network.facets)
        for (int side = 0; side < 2; ++side)
            checkTrue(network.cellPores[f.cells[side]] == static_cast<long>(f.pores[side]), "random packing: facet tetrahedra in their pores");

    // the grid of the network data and the same grid read back from a DGF file give the same flow
    Dumux::PoreNetwork::GridManager<3> dataGridManager, dgfGridManager;
    dataGridManager.init(networkData(network));
    Dumux::PoreNetwork::writeDgf("spherepacking.dgf", dataGridManager.grid().leafGridView(), *dataGridManager.getGridData());
    dgfGridManager.init("Dgf");
    const auto flow = solveFlow(dataGridManager, network);
    check(flow.outflow, flow.inflow, 1e-10, "random packing: outflow equals inflow");
    check(flow.inflow, flow.networkFlux, 1e-10, "random packing: model flux equals the flux from the transmissibilities");
    check(solveFlow(dgfGridManager, network).inflow, flow.inflow, 1e-14, "random packing: same flux for the grid read from a DGF file");

    // given transmissibilities, the Chareyre law where none is given
    using Given = Dumux::Properties::TTag::SpherePackingGivenTransmissibility;
    check(solveFlow<Given>(network).inflow, flow.inflow, 1e-12, "given transmissibility: none given, Chareyre law");
    auto doubled = network;
    for (auto& throat : doubled.throats)
        throat.transmissibility = 2.0*throatTransmissibility(throat);
    check(solveFlow<Given>(doubled).inflow, 2.0*flow.inflow, 1e-12, "given transmissibility: doubled transmissibilities double the flux");
    std::mt19937 gen(11);
    std::uniform_real_distribution<double> logFactor(std::log(0.25), std::log(4.0));
    auto given = network;
    for (auto& throat : given.throats)
        throat.transmissibility = std::exp(logFactor(gen))*throatTransmissibility(throat);
    const auto givenFlow = solveFlow<Given>(given);
    check(givenFlow.inflow, givenFlow.networkFlux, 1e-10, "given transmissibility: model flux equals the flux from the transmissibilities");
}

} // end anonymous namespace

int main(int argc, char** argv)
{
    Dumux::initialize(argc, argv);
    Dumux::Parameters::init([](auto& tree) {
        tree["Dgf.Grid.File"] = "spherepacking.dgf";
        tree["Grid.PoreGeometry"] = "Cube";
        tree["Problem.EnableGravity"] = "false";
        tree["Problem.InletPressure"] = "1";
        tree["Problem.OutletPressure"] = "0";
        tree["Component.LiquidDensity"] = "1000";
        tree["Component.LiquidDynamicViscosity"] = "1e-3";
    });

    testCubicLattice();
    testRandomPacking();

    if (numFailures > 0)
        DUNE_THROW(Dune::Exception, numFailures << " checks failed");
    std::cout << "all checks passed\n";
    return 0;
}
