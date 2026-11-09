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
#include <dumux/io/grid/porenetwork/dgfwriter.hh>
#include <dumux/porenetwork/common/boundaryflux.hh>
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

enum class GridSource { dgfFile, networkData };

//! DGF file of the network, written by the DGF writer of the pore-network grids
void writeNetworkDgf(const std::string& fileName, const Network<double>& network)
{
    Dumux::PoreNetwork::GridManager<3> gridManager;
    gridManager.init(networkData(network));
    Dumux::PoreNetwork::writeDgf(fileName, gridManager.grid().leafGridView(), *gridManager.getGridData());
}

template<class TypeTag = Dumux::Properties::TTag::SpherePackingPermeability>
FlowResult solveFlow(const std::string& paramGroup, const Network<double>& network, int axis,
                     GridSource source = GridSource::dgfFile)
{
    using namespace Dumux;

    PoreNetwork::GridManager<3> gridManager;
    if (source == GridSource::networkData)
        gridManager.init(networkData(network), paramGroup);
    else
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
        const double transmissibility = throatTransmissibility(throat);
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
    writeNetworkDgf("lattice.dgf", network);

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
    writeNetworkDgf("random.dgf", network);

    const auto flow = solveFlow("Random", network, 0);
    check(flow.modelOutflux, flow.modelInflux, 1e-10, "random packing: outflow equals inflow");
    check(flow.modelInflux, flow.networkFlux, 1e-10, "random packing: model flux equals network flux");

    const auto fromData = solveFlow("Random", network, 0, GridSource::networkData);
    check(fromData.modelInflux, flow.modelInflux, 1e-12, "random packing: grid from network data, same flux as from the DGF file");
    double maxDifference = 0.0;
    for (std::size_t i = 0; i < network.pores.size(); ++i)
        maxDifference = std::max(maxDifference, std::abs(fromData.pressure[i] - flow.pressure[i]));
    check(1.0 + maxDifference, 1.0, 1e-12, "random packing: grid from network data, same pressures as from the DGF file");

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

void testGivenTransmissibility()
{
    using GivenTransmissibility = Dumux::Properties::TTag::SpherePackingGivenTransmissibility;
    const auto [centers, radii, boxLength] = randomPacking(1200, 0.25, 7);
    ExtractionRegion<double> region;
    region.lower = 1.5e-3;
    region.upper = boxLength - 1.5e-3;
    auto network = extractNetwork(regularTriangulation(centers, radii), centers, radii, region);
    // inlet and outlet only
    for (auto& pore : network.pores)
        pore.label = pore.label == 1 || pore.label == 2 ? pore.label : -1;

    const auto chareyre = solveFlow("Given", network, 0, GridSource::networkData);
    const auto fromGeometry = solveFlow<GivenTransmissibility>("Given", network, 0, GridSource::networkData);
    check(fromGeometry.modelInflux, chareyre.modelInflux, 1e-12, "given transmissibility: none given, Chareyre law");

    std::mt19937 gen(11);
    std::uniform_real_distribution<double> logFactor(std::log(0.25), std::log(4.0));
    auto given = network;
    for (auto& throat : given.throats)
        throat.transmissibility = std::exp(logFactor(gen))*throatTransmissibility(throat);
    const auto flow = solveFlow<GivenTransmissibility>("Given", given, 0, GridSource::networkData);
    check(flow.modelOutflux, flow.modelInflux, 1e-10, "given transmissibility: outflow equals inflow");
    check(flow.modelInflux, flow.networkFlux, 1e-10, "given transmissibility: model flux equals network flux");
    checkTrue(std::abs(flow.modelInflux/chareyre.modelInflux - 1.0) > 1e-3, "given transmissibility: flux differs from Chareyre");

    auto doubled = network;
    for (auto& throat : doubled.throats)
        throat.transmissibility = 2.0*throatTransmissibility(throat);
    const auto twice = solveFlow<GivenTransmissibility>("Given", doubled, 0, GridSource::networkData);
    check(twice.modelInflux, 2.0*chareyre.modelInflux, 1e-12, "given transmissibility: doubled transmissibilities double the flux");
}

} // end anonymous namespace

int main(int argc, char** argv)
{
    Dumux::initialize(argc, argv);
    Dumux::Parameters::init([](auto& tree) {
        tree["Lattice.Grid.File"] = "lattice.dgf";
        tree["Random.Grid.File"] = "random.dgf";
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
    testGivenTransmissibility();

    if (numFailures > 0)
        DUNE_THROW(Dune::Exception, numFailures << " checks failed");
    std::cout << "all checks passed\n";
    return 0;
}
