// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Dimensionless Boussinesq dissolution test — pressure-based formulation.
 *
 * One-sided Rayleigh-Bénard convection using the standard OnePNC model with
 * BoussinesqCVFEDarcyLaw.  Pressure reference is fixed weakly by a Robin
 * penalty BC at the top (Nitsche trick); replaceCompEqIdx = 0 makes eq 0 the
 * total mass balance (divergence-free velocity constraint).
 *
 * Mirrors main_vorticity.cc (perturbed IC, n_hat/tip_depth diagnostics) so the
 * same experiment (same Ra, PerturbAmplitude, PerturbSeed) can be run through
 * both formulations and cross-checked against each other.
 */

#include <config.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/fvassembler.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#include <dumux/discretization/evalgradients.hh>
#include <dumux/discretization/cvfe/elementsolution.hh>

#include "properties_pressure.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTag = Properties::TTag::BoussinesqPressureTest;

    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    Parameters::init(argc, argv, [](auto&){}, "params.input");
    GetPropType<TypeTag, Properties::FluidSystem>::init();  // reads Ra

    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    const auto& leafGridView = gridManager.grid().leafGridView();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    using Scalar  = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    static constexpr int dimWorld = GridGeometry::GridView::dimensionworld;

    // primary variable index for solute concentration (mass fraction of component 1)
    static constexpr int concentrationIdx = FluidSystem::soluteIdx;

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    problem->applyInitialSolution(x);

    const Scalar Ra = getParam<Scalar>("DimensionlessNumbers.Ra");

    // Optionally start from the analytical diffusive profile
    const Scalar tcrit  = 146.0 / Ra;
    const Scalar tStart = getParam<Scalar>("TimeLoop.TStart", tcrit);
    const Scalar perturbAmplitude = getParam<Scalar>("Problem.PerturbAmplitude", 0.0);
    const int    perturbSeed      = getParam<int>("Problem.PerturbSeed", 1);
    if (tStart > 0.0)
    {
        const Scalar zMax = gridGeometry->bBoxMax()[dimWorld-1];

        // One N(0,1) draw per x-coordinate, memoized so every vertex in the
        // same z-column (same x, structured grid) reuses the identical draw.
        // Same recipe as main_vorticity.cc.
        std::mt19937 rng(static_cast<unsigned int>(perturbSeed));
        std::normal_distribution<Scalar> noiseDist(0.0, 1.0);
        std::map<Scalar, Scalar> noiseByX;
        const auto noiseAt = [&](Scalar xCoord) -> Scalar {
            auto it = noiseByX.find(xCoord);
            if (it == noiseByX.end())
                it = noiseByX.emplace(xCoord, noiseDist(rng)).first;
            return it->second;
        };

        std::size_t numClipped = 0;
        std::size_t numVertices = 0;
        for (const auto& vertex : vertices(leafGridView))
        {
            const auto idx = gridGeometry->dofMapper().index(vertex);
            const auto pos = vertex.geometry().center();
            const Scalar z   = pos[dimWorld-1];
            const Scalar eta = (zMax - z) / (2.0 * std::sqrt(tStart / Ra));
            Scalar C = std::erfc(eta);

            if (perturbAmplitude > 0.0)
            {
                const Scalar envelope = std::exp(-eta*eta);
                C += perturbAmplitude * envelope * noiseAt(pos[0]);
                ++numVertices;
                if (C < 0.0 || C > 1.0)
                    ++numClipped;
                C = std::min(std::max(C, 0.0), 1.0);
            }
            x[idx][concentrationIdx] = C;
        }

        if (perturbAmplitude > 0.0 && mpiHelper.rank() == 0)
            std::cout << "Perturbed IC: amplitude = " << perturbAmplitude
                      << ", seed = " << perturbSeed
                      << ", clipped " << numClipped << "/" << numVertices
                      << " vertices to [0,1]\n";
    }

    auto xOld = x;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // --- Sherwood number: dimensionless dissolution flux at the top boundary ---
    const Scalar zMax = gridGeometry->bBoxMax()[dimWorld-1];

    Scalar domainTopArea = 1.0;
    for (int d = 0; d < dimWorld-1; ++d)
        domainTopArea *= gridGeometry->bBoxMax()[d] - gridGeometry->bBoxMin()[d];

    auto computeSherwood = [&]() -> Scalar {
        Scalar totalFlux = 0.0;
        for (const auto& element : elements(leafGridView, Dune::Partitions::interior))
        {
            auto fvGeometry = localView(*gridGeometry);
            fvGeometry.bindElement(element);

            bool hasTopFace = false;
            for (const auto& scvf : scvfs(fvGeometry))
                if (scvf.boundary() && scvf.center()[dimWorld-1] > zMax - 1e-8)
                    { hasTopFace = true; break; }
            if (!hasTopFace) continue;

            const auto elemSol = elementSolution(element, x, *gridGeometry);
            const auto geo = element.geometry();

            for (const auto& scvf : scvfs(fvGeometry))
            {
                if (!scvf.boundary() || scvf.center()[dimWorld-1] < zMax - 1e-8)
                    continue;

                const auto grads = evalGradients(element, geo, *gridGeometry, elemSol, scvf.center());
                totalFlux += (grads[concentrationIdx] * scvf.unitOuterNormal()) * scvf.area();
            }
        }
        totalFlux = leafGridView.comm().sum(totalFlux);
        return totalFlux / domainTopArea;
    };

    auto computeTotalMass = [&]() -> Scalar {
        Scalar mass = 0.0;
        for (const auto& element : elements(leafGridView, Dune::Partitions::interior))
        {
            auto fvGeometry = localView(*gridGeometry);
            fvGeometry.bindElement(element);
            for (const auto& scv : scvs(fvGeometry))
                mass += x[scv.dofIndex()][concentrationIdx] * scv.volume();
        }
        return leafGridView.comm().sum(mass);
    };

    // --- Dominant finger count n_hat(t) and finger tip depth: see
    // main_vorticity.cc for the derivation/references (Riaz et al. 2006 eq.
    // 3.10/3.11). Identical formulas -- concentration is
    // a primary variable in both formulations, so these diagnostics don't
    // depend on which one computed it. ---
    constexpr Scalar pi = 3.14159265358979323846;
    const Scalar xMin = gridGeometry->bBoxMin()[0];
    const Scalar xMax = gridGeometry->bBoxMax()[0];
    const Scalar Lx   = xMax - xMin;
    int nx = 0;
    for (const auto c : getParam<std::vector<int>>("Grid.Cells0"))
        nx += c;
    const Scalar dxCol = Lx / nx;

    // Returns {n_hat, totalVorticityEnergy}; see the vorticity-formulation
    // main.cc for why the energy needs to be reported alongside n_hat.
    auto computeDominantWavenumberAndEnergy = [&]() -> std::pair<Scalar, Scalar> {
        std::vector<Scalar> g(nx, 0.0);
        for (const auto& element : elements(leafGridView, Dune::Partitions::interior))
        {
            auto fvGeometry = localView(*gridGeometry);
            fvGeometry.bindElement(element);
            const auto elemSol = elementSolution(element, x, *gridGeometry);
            const auto geo = element.geometry();
            const auto center = geo.center();

            const auto grads = evalGradients(element, geo, *gridGeometry, elemSol, center);
            const Scalar omega = -grads[concentrationIdx][0]; // omega = -dC/dx

            int col = static_cast<int>(std::floor((center[0] - xMin) / dxCol));
            col = std::min(std::max(col, 0), nx - 1);
            g[col] += omega * geo.volume();
        }
        for (auto& v : g)
            v = leafGridView.comm().sum(v) / dxCol;

        const int nModes = nx / 2;
        Scalar numerator = 0.0, denominator = 0.0;
        for (int n = 1; n <= nModes; ++n)
        {
            const Scalar k = 2.0*pi*n / Lx;
            Scalar reSum = 0.0, imSum = 0.0;
            for (int i = 0; i < nx; ++i)
            {
                const Scalar xi = xMin + (i + 0.5) * dxCol;
                reSum += g[i] * std::cos(k * xi) * dxCol;
                imSum -= g[i] * std::sin(k * xi) * dxCol;
            }
            const Scalar E = reSum*reSum + imSum*imSum;
            numerator   += k * E;
            denominator += E;
        }
        const Scalar nHat = denominator > 0.0 ? (numerator / denominator) / (2.0*pi) : 0.0;
        return {nHat, denominator};
    };

    const Scalar tipThreshold = getParam<Scalar>("Problem.TipThreshold", 0.01);
    auto computeTipDepth = [&]() -> Scalar {
        Scalar minZ = zMax;
        for (const auto& vertex : vertices(leafGridView))
        {
            const auto idx = gridGeometry->dofMapper().index(vertex);
            if (x[idx][concentrationIdx] >= tipThreshold)
                minZ = std::min(minZ, vertex.geometry().center()[dimWorld-1]);
        }
        minZ = leafGridView.comm().min(minZ);
        return zMax - minZ;
    };

    std::ofstream shFile(problem->name() + "_sherwood.csv");
    shFile << std::scientific << std::setprecision(6);
    shFile << "time,Sh_grad,F_grad,Sh_mass,F_mass\n";

    std::ofstream fingerFile(problem->name() + "_finger.csv");
    fingerFile << std::scientific << std::setprecision(6);
    fingerFile << "time,n_hat,vorticity_energy,tip_depth\n";

    const auto tEnd         = getParam<Scalar>("TimeLoop.TEnd");
    const auto dtInit       = getParam<Scalar>("TimeLoop.DtInitial");
    const auto maxDt        = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    const auto vtkInterval  = getParam<Scalar>("Output.VtkOutputInterval", 10.0);

    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    IOFields::initOutputModule(vtkWriter);
    vtkWriter.write(0.0);

    const auto writeSherwood = [&](Scalar t, Scalar massFluxPerWidth = 0.0) {
        const Scalar Sh_grad = computeSherwood();
        const Scalar Sh_mass = massFluxPerWidth * Ra;
        shFile << t << "," << Sh_grad << "," << Sh_grad/Ra
                    << "," << Sh_mass << "," << massFluxPerWidth << "\n";
        shFile.flush();
    };

    const auto writeFinger = [&](Scalar t) {
        const auto [nHat, vorticityEnergy] = computeDominantWavenumberAndEnergy();
        fingerFile << t << "," << nHat << "," << vorticityEnergy << "," << computeTipDepth() << "\n";
        fingerFile.flush();
    };

    Scalar massOld = computeTotalMass();
    writeSherwood(tStart);
    writeFinger(tStart);

    Scalar nextVtkTime = tStart + vtkInterval;
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(tStart, dtInit, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);

    using LinearSolver = Dumux::AMGBiCGSTABIstlSolver<SeqLinearSolverTraits,
                       LinearAlgebraTraitsFromAssembler<Assembler>>;

    auto linearSolver = std::make_shared<LinearSolver>();

    NewtonSolver<Assembler, LinearSolver> nonLinearSolver(assembler, linearSolver);

    timeLoop->start(); do
    {
        nonLinearSolver.solve(x, *timeLoop);

        const Scalar dt = timeLoop->timeStepSize();
        const Scalar massNew = computeTotalMass();
        const Scalar massFluxPerWidth = (massNew - massOld) / (dt * domainTopArea);
        massOld = massNew;

        xOld = x;
        gridVariables->advanceTimeStep();
        timeLoop->advanceTimeStep();

        if (timeLoop->time() >= nextVtkTime - 1e-8*vtkInterval || timeLoop->finished())
        {
            vtkWriter.write(timeLoop->time());
            nextVtkTime += vtkInterval;
        }

        writeSherwood(timeLoop->time(), massFluxPerWidth);
        writeFinger(timeLoop->time());

        timeLoop->reportTimeStep();
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/false);

    return 0;
}
