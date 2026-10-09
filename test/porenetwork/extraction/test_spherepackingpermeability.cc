// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkModels
 * \brief Test of single-phase flow through the pore network of a sphere packing
 *
 * Simple cubic lattice: the permeability of the network is A R_h^2/(2 a^2) with the fluid area A
 * and hydraulic radius R_h of a face and the lattice constant a. Random packing: inflow equals
 * outflow, and the flux computed by the model equals the flux from the extracted throat parameters.
 */
#include <config.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <optional>
#include <numbers>
#include <random>
#include <string>
#include <tuple>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/io/grid/porenetwork/gridmanager.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/porenetwork/common/boundaryflux.hh>
#include <dumux/porenetwork/extraction/porescalefinitevolume.hh>
#include <dumux/porenetwork/extraction/spherepackingtriangulation.hh>

#include "properties.hh"

namespace {

using namespace Dumux::PoreNetwork::SpherePacking;
using P = Point<double>;
constexpr double pi = std::numbers::pi;

int numFailures = 0;

void check(double value, double reference, double tolerance, const std::string& name)
{
    const double error = std::abs(value - reference)/std::max(std::abs(reference), 1e-300);
    std::cout << name << ": " << value << " vs " << reference << " (relative error " << error << ")\n";
    if (!(error <= tolerance))
    {
        std::cout << "FAILED " << name << "\n";
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

struct FlowResult
{
    double modelInflux; //!< volume flux into the network through the inlet pores computed by the model
    double modelOutflux; //!< volume flux out of the network through the outlet pores computed by the model
    double networkFlux; //!< volume flux through the inlet pores from the throat parameters and the pressures
    double inletPosition, outletPosition; //!< mean coordinate of the inlet and outlet pores along the flow axis
    std::vector<double> pressure; //!< pressure of each pore of the network
};

FlowResult solveFlow(const std::string& paramGroup, const Network<double>& network, int axis)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::SpherePackingPermeability;

    PoreNetwork::GridManager<3> gridManager;
    gridManager.init(paramGroup);
    const auto& gridView = gridManager.grid().leafGridView();
    const auto gridData = gridManager.getGridData();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(gridView, *gridData);
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;
    auto spatialParams = std::make_shared<SpatialParams>(gridGeometry, *gridData);
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
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    LinearSolver linearSolver;
    linearSolver.solve(*A, x, *r);
    gridVariables->update(x);

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    const double density = FluidSystem::density(0.0, 0.0);
    const double viscosity = FluidSystem::viscosity(0.0, 0.0);
    const auto boundaryFlux = PoreNetwork::BoundaryFlux(*gridVariables, assembler->localResidual(), x);

    FlowResult result;
    // the boundary flux is positive for flow out of the network
    result.modelInflux = -boundaryFlux.getFlux(std::vector<int>{problem->inletLabel()}).totalFlux[0]/density;
    result.modelOutflux = boundaryFlux.getFlux(std::vector<int>{problem->outletLabel()}).totalFlux[0]/density;

    // flux from the network data, matching the grid vertices to the pores by position
    std::map<std::array<long, 3>, std::size_t> poreAt;
    const auto key = [](const P& p) {
        return std::array<long, 3>{std::lround(p[0]*1e12), std::lround(p[1]*1e12), std::lround(p[2]*1e12)};
    };
    for (std::size_t i = 0; i < network.pores.size(); ++i)
        poreAt[key(network.pores[i].position)] = i;
    std::vector<double> pressure(network.pores.size());
    for (const auto& vertex : vertices(gridView))
        pressure[poreAt.at(key(vertex.geometry().center()))] = x[gridGeometry->vertexMapper().index(vertex)][0];

    result.pressure = pressure;
    result.networkFlux = 0.0;
    for (const auto& throat : network.throats)
    {
        const auto i = throat.pores[0], j = throat.pores[1];
        const bool inletI = network.pores[i].label == problem->inletLabel();
        const bool inletJ = network.pores[j].label == problem->inletLabel();
        if (inletI == inletJ)
            continue;
        const double transmissibility = PoreNetwork::TransmissibilityChareyre<double>::singlePhaseTransmissibility(
            throat.fluidArea, throat.hydraulicRadius, throat.length);
        const double dp = inletI ? pressure[i] - pressure[j] : pressure[j] - pressure[i];
        result.networkFlux += transmissibility/viscosity*dp;
    }

    double inletSum = 0.0, outletSum = 0.0;
    int numInlet = 0, numOutlet = 0;
    for (const auto& pore : network.pores)
    {
        if (pore.label == problem->inletLabel()) { inletSum += pore.position[axis]; ++numInlet; }
        if (pore.label == problem->outletLabel()) { outletSum += pore.position[axis]; ++numOutlet; }
    }
    result.inletPosition = inletSum/numInlet;
    result.outletPosition = outletSum/numOutlet;
    return result;
}

void testCubicLattice()
{
    const int n = 6;
    const double a = 1e-3, r = 0.45e-3;
    std::vector<P> centers;
    std::vector<double> radii;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
            {
                centers.push_back(P{i*a, j*a, k*a});
                radii.push_back(r);
            }

    ExtractionRegion<double> region;
    region.lower = -0.1*a;
    region.upper = (n - 1 + 0.1)*a;
    region.mergeDistance = 1e-9*a;
    auto network = extractNetwork(regularTriangulation(centers, radii), centers, radii, region);
    for (auto& pore : network.pores)
        pore.label = pore.position[0] < a ? 1 : (pore.position[0] > (n - 2)*a ? 2 : -1);
    writeDgf("lattice.dgf", network);

    const auto flow = solveFlow("Lattice", network, 0);
    const double viscosity = Dumux::getParam<double>("Component.LiquidDynamicViscosity");
    const double dp = Dumux::getParam<double>("Problem.InletPressure") - Dumux::getParam<double>("Problem.OutletPressure");
    const double length = flow.outletPosition - flow.inletPosition;
    const double area = (n - 1)*a*(n - 1)*a;
    const double permeability = viscosity*flow.modelInflux*length/(area*dp);

    const double fluidArea = a*a - pi*r*r;
    const double hydraulicRadius = (a*a*a/3.0 - 4.0*pi*r*r*r/9.0)/(4.0*pi*r*r/3.0);
    check(length, (n - 2)*a, 1e-12, "lattice: distance of inlet and outlet pores");
    check(permeability, fluidArea*hydraulicRadius*hydraulicRadius/(2.0*a*a), 1e-10, "lattice: permeability");
    check(flow.modelOutflux, flow.modelInflux, 1e-10, "lattice: outflow equals inflow");
}

void testRandomPacking()
{
    // random sequential addition of non-overlapping spheres, largest first
    const int n = 2000;
    const double solidFraction = 0.25;
    std::mt19937 gen(3);
    std::uniform_real_distribution<double> radius(0.4e-3, 0.6e-3), unit(0.0, 1.0);
    std::vector<double> radii(n);
    for (auto& r : radii)
        r = radius(gen);
    std::sort(radii.begin(), radii.end(), std::greater<double>{});
    double solid = 0.0;
    for (const double r : radii)
        solid += 4.0/3.0*pi*r*r*r;
    const double boxLength = std::cbrt(solid/solidFraction);
    std::vector<P> centers;
    for (const double r : radii)
    {
        while (true)
        {
            const P x{r + (boxLength - 2*r)*unit(gen), r + (boxLength - 2*r)*unit(gen), r + (boxLength - 2*r)*unit(gen)};
            if (std::all_of(centers.begin(), centers.end(), [&, i = std::size_t(0)](const P& c) mutable {
                    return (x - c).two_norm() >= r + radii[i++]; }))
            {
                centers.push_back(x);
                break;
            }
        }
    }

    ExtractionRegion<double> region;
    region.lower = 2e-3;
    region.upper = boxLength - 2e-3;
    const auto network = extractNetwork(regularTriangulation(centers, radii), centers, radii, region);
    writeDgf("random.dgf", network);

    const auto flow = solveFlow("Random", network, 0);
    check(flow.modelOutflux, flow.modelInflux, 1e-10, "random packing: outflow equals inflow");
    check(flow.modelInflux, flow.networkFlux, 1e-10, "random packing: model flux equals network flux");

    const double viscosity = Dumux::getParam<double>("Component.LiquidDynamicViscosity");
    const double dp = Dumux::getParam<double>("Problem.InletPressure") - Dumux::getParam<double>("Problem.OutletPressure");
    const double extent = region.upper[1] - region.lower[1];
    std::cout << "random packing: permeability "
              << viscosity*flow.modelInflux*(flow.outletPosition - flow.inletPosition)/(extent*extent*dp) << " m^2\n";
}

// random sequential addition of non-overlapping spheres in a box, largest first
std::tuple<std::vector<P>, std::vector<double>, double> randomPacking(int n, double solidFraction, unsigned int seed)
{
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> radius(0.4e-3, 0.6e-3), unit(0.0, 1.0);
    std::vector<double> radii(n);
    for (auto& r : radii)
        r = radius(gen);
    std::sort(radii.begin(), radii.end(), std::greater<double>{});
    double solid = 0.0;
    for (const double r : radii)
        solid += 4.0/3.0*pi*r*r*r;
    const double boxLength = std::cbrt(solid/solidFraction);
    std::vector<P> centers;
    for (const double r : radii)
    {
        while (true)
        {
            const P x{r + (boxLength - 2*r)*unit(gen), r + (boxLength - 2*r)*unit(gen), r + (boxLength - 2*r)*unit(gen)};
            if (std::all_of(centers.begin(), centers.end(), [&, i = std::size_t(0)](const P& c) mutable {
                    return (x - c).two_norm() >= r + radii[i++]; }))
            {
                centers.push_back(x);
                break;
            }
        }
    }
    return {centers, radii, boxLength};
}

void testWalls()
{
    const auto [centers, radii, boxLength] = randomPacking(1500, 0.25, 5);
    Walls<double> walls;
    walls.lower = 0.0;
    walls.upper = boxLength;
    walls.numSpheres = centers.size();
    const auto network = extractNetwork(regularTriangulation(centers, radii, walls), centers, radii, walls);

    // the pore-network model with the pores touching the walls x = min and x = max as inlet and outlet
    auto labelled = network;
    for (auto& pore : labelled.pores)
        pore.label = pore.label > 0 && (pore.label & 1) ? 1 : (pore.label > 0 && (pore.label & 2) ? 2 : -1);
    writeDgf("walls.dgf", labelled);
    const auto flow = solveFlow("Walls", labelled, 0);
    check(flow.modelOutflux, flow.modelInflux, 1e-10, "walls: outflow equals inflow");

    // the pore-scale flow solver gives the same pressures
    const double viscosity = Dumux::getParam<double>("Component.LiquidDynamicViscosity");
    const double dp = Dumux::getParam<double>("Problem.InletPressure") - Dumux::getParam<double>("Problem.OutletPressure");
    PoreScaleFlow<double> poreScaleFlow(network, viscosity, {Dumux::getParam<double>("Problem.InletPressure"),
                                                              Dumux::getParam<double>("Problem.OutletPressure"),
                                                              std::nullopt, std::nullopt, std::nullopt, std::nullopt});
    const auto& pressure = poreScaleFlow.solve(std::vector<double>(network.pores.size(), 0.0));
    double maxDifference = 0.0;
    for (std::size_t i = 0; i < pressure.size(); ++i)
        maxDifference = std::max(maxDifference, std::abs(pressure[i] - flow.pressure[i]));
    check(1.0 + maxDifference/dp, 1.0, 1e-10, "walls: pore-scale flow pressures equal the pore-network model");

    // uniform compression drained at x = min and x = max: p = mu rate x (L - x)/(2 k) in the mean
    const double permeability = viscosity*flow.modelInflux*boxLength/(boxLength*boxLength*dp);
    const double rate = 1e-3;
    std::vector<double> volumeRates;
    for (const auto& pore : network.pores)
        volumeRates.push_back(-rate*pore.bulkVolume);
    PoreScaleFlow<double> drained(network, viscosity, {0.0, 0.0, std::nullopt, std::nullopt, std::nullopt, std::nullopt});
    const auto& p = drained.solve(volumeRates);
    double middle = 0.0;
    int numMiddle = 0;
    for (std::size_t i = 0; i < p.size(); ++i)
        if (std::abs(network.pores[i].position[0] - 0.5*boxLength) < 0.05*boxLength)
        {
            middle += p[i];
            ++numMiddle;
        }
    check(middle/numMiddle, viscosity*rate*boxLength*boxLength/(8.0*permeability), 0.1,
          "walls: mean pressure in the middle of a uniformly compressed sample");
    checkTrue(*std::min_element(p.begin(), p.end()) >= 0.0, "walls: no suction under compression");
}

// three steps of uniform compression with a new triangulation before the third
void testPoreScaleFiniteVolume()
{
    const auto [centers, radii, boxLength] = randomPacking(800, 0.25, 6);
    Walls<double> walls;
    walls.lower = 0.0;
    walls.upper = boxLength;
    walls.numSpheres = centers.size();
    const double viscosity = 1e-3, rate = 1.0, dt = 1e-5;
    const std::array<std::optional<double>, 6> wallPressure{0.0, 0.0, std::nullopt, std::nullopt, std::nullopt, std::nullopt};
    const P center(0.5*boxLength);
    const auto compress = [&](auto x, auto w, int steps) {
        const double factor = std::pow(1.0 - rate*dt, steps);
        for (auto& c : x)
            c = center + factor*(c - center);
        w.lower = center + factor*(w.lower - center);
        w.upper = center + factor*(w.upper - center);
        return std::make_pair(x, w);
    };

    PoreScaleFiniteVolume<double> fluid(viscosity, wallPressure);
    fluid.remesh(centers, radii, walls);
    for (int step = 1; step <= 3; ++step)
    {
        if (step == 3)
        {
            const auto [x, w] = compress(centers, walls, 2);
            fluid.remesh(x, radii, w);
        }
        const auto [x, w] = compress(centers, walls, step);
        const auto [xOld, wOld] = compress(centers, walls, step - 1);
        fluid.update(x, w, dt);

        // the same step from scratch: triangulation at the previous positions, one update; in the second step
        // the network is that of the initial positions, which differs from the new one by the strain rate dt
        PoreScaleFiniteVolume<double> reference(viscosity, wallPressure);
        reference.remesh(xOld, radii, wOld);
        reference.update(x, w, dt);
        double difference = 0.0, scale = 0.0;
        for (std::size_t i = 0; i < reference.pressure().size(); ++i)
        {
            difference = std::max(difference, std::abs(fluid.pressure()[i] - reference.pressure()[i]));
            scale = std::max(scale, std::abs(reference.pressure()[i]));
        }
        check(1.0 + difference/scale, 1.0, step == 2 ? 1e-3 : 1e-10, "pore-scale finite volume: pressures of step " + std::to_string(step));

        const auto& forces = fluid.forces();
        checkTrue(forces.size() == centers.size() + 6, "pore-scale finite volume: forces on spheres and walls");
        checkTrue(forces[centers.size()].two_norm() == 0.0 && forces[centers.size() + 1].two_norm() == 0.0,
                  "pore-scale finite volume: no force on the drained walls");
    }
}

// compressible fluid: storage, volume sources, pressure transfer to a new network, incompressible limit
void testCompressible()
{
    const auto [centers, radii, boxLength] = randomPacking(600, 0.25, 7);
    Walls<double> walls;
    walls.lower = 0.0;
    walls.upper = boxLength;
    walls.numSpheres = centers.size();
    const double viscosity = 1e-3, dt = 1e-5, bulkModulus = 1e8, source = 2e-3;
    const std::array<std::optional<double>, 6> closed{};

    // a closed box with a uniform volume source per void volume: dp/dt = K q in every pore
    PoreScaleFiniteVolume<double> fluid(viscosity, closed, bulkModulus, dt);
    fluid.remesh(centers, radii, walls);
    for (int step = 1; step <= 3; ++step)
    {
        if (step == 3)
            fluid.remesh(centers, radii, walls);
        fluid.update(centers, walls, dt, source);
        const auto [low, high] = std::minmax_element(fluid.pressure().begin(), fluid.pressure().end());
        check(*low, step*bulkModulus*source*dt, 1e-8, "compressible: uniform pressure from a volume source (min), step " + std::to_string(step));
        check(*high, step*bulkModulus*source*dt, 1e-8, "compressible: uniform pressure from a volume source (max), step " + std::to_string(step));
    }
    const auto& forces = fluid.forces();
    double sphereForce = 0.0;
    for (std::size_t i = 0; i < centers.size(); ++i)
        sphereForce = std::max(sphereForce, forces[i].two_norm());
    checkTrue(sphereForce < 1e-12*3*bulkModulus*source*dt*boxLength*boxLength, "compressible: no force on the spheres for a uniform pressure");

    // a very stiff fluid behaves as the incompressible one
    const std::array<std::optional<double>, 6> drained{0.0, 0.0, std::nullopt, std::nullopt, std::nullopt, std::nullopt};
    const P center(0.5*boxLength);
    auto moved = centers;
    for (auto& x : moved)
        x = center + (1.0 - dt)*(x - center);
    auto movedWalls = walls;
    movedWalls.lower = center + (1.0 - dt)*(walls.lower - center);
    movedWalls.upper = center + (1.0 - dt)*(walls.upper - center);
    PoreScaleFiniteVolume<double> incompressible(viscosity, drained);
    PoreScaleFiniteVolume<double> stiff(viscosity, drained, 1e20, dt);
    incompressible.remesh(centers, radii, walls);
    stiff.remesh(centers, radii, walls);
    incompressible.update(moved, movedWalls, dt);
    stiff.update(moved, movedWalls, dt);
    double difference = 0.0, scale = 0.0;
    for (std::size_t i = 0; i < stiff.pressure().size(); ++i)
    {
        difference = std::max(difference, std::abs(stiff.pressure()[i] - incompressible.pressure()[i]));
        scale = std::max(scale, std::abs(incompressible.pressure()[i]));
    }
    check(1.0 + difference/scale, 1.0, 1e-8, "compressible: incompressible limit");
}

} // end anonymous namespace

int main(int argc, char** argv)
{
    Dumux::initialize(argc, argv);
    Dumux::Parameters::init([](auto& tree) {
        tree["Lattice.Grid.File"] = "lattice.dgf";
        tree["Random.Grid.File"] = "random.dgf";
        tree["Walls.Grid.File"] = "walls.dgf";
        tree["Grid.PoreGeometry"] = "Cube";
        tree["Grid.Sanitize"] = "false";
        tree["Problem.Name"] = "spherepacking";
        tree["Problem.EnableGravity"] = "false";
        tree["Problem.InletLabel"] = "1";
        tree["Problem.OutletLabel"] = "2";
        tree["Problem.InletPressure"] = "1";
        tree["Problem.OutletPressure"] = "0";
        tree["Component.LiquidDensity"] = "1000";
        tree["Component.LiquidDynamicViscosity"] = "1e-3";
    });

    testCubicLattice();
    testRandomPacking();
    testWalls();
    testPoreScaleFiniteVolume();
    testCompressible();

    if (numFailures > 0)
        DUNE_THROW(Dune::Exception, numFailures << " checks failed");
    std::cout << "all checks passed\n";
    return 0;
}
