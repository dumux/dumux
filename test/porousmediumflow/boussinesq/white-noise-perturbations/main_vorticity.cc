// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Dimensionless Boussinesq dissolution test: one-sided Rayleigh-Bénard convection,
 *        with a (localized, white-noise) perturbed initial condition.
 *
 * With Problem.PerturbAmplitude left at its default (0), this reduces exactly to the
 * unperturbed erfc base state used for the dissolution-flux (Sherwood) comparison
 * against the pressure-based formulation.
 *
 * With Problem.PerturbAmplitude > 0, a localized white-noise perturbation is
 * superimposed on the erfc base state at t = tStart (one N(0,1) draw per x,
 * shared down the z-column, scaled by exp(-eta^2) so it is confined to the
 * diffusive boundary layer -- see Riaz, Hesse, Tchelepi & Orr (2006) sec. 3.2).
 *
 * Two diagnostics are written to <name>_finger.csv at every accepted timestep:
 *   - n_hat(t): dominant finger count, the energy-weighted mean wavenumber of
 *     the depth-integrated vorticity field omega = -dC/dx (their eq. 3.10/3.11),
 *     converted from radian wavenumber to a mode/finger count via /(2*pi)
 *     (valid since the domain width here is 1, matching their aspect ratio A=1).
 *   - tip_depth(t): zMax minus the shallowest z at which C >= Problem.TipThreshold
 *     for any x (i.e. how far the most advanced finger has descended from the
 *     top wall); the threshold is not reported by the paper and must be chosen
 *     and documented explicitly (default 0.01).
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

#include "properties_vorticity.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTag = Properties::TTag::TYPETAG;

    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    constexpr int dimW = GetPropType<TypeTag, Properties::Grid>::dimensionworld;
    const std::string defaultParamFile = (dimW == 3) ? "params3d.input" : "params.input";
    Parameters::init(argc, argv, [](auto&){}, defaultParamFile);
    GetPropType<TypeTag, Properties::FluidSystem>::init();

    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    const auto& leafGridView = gridManager.grid().leafGridView();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    static constexpr int dimWorld = GridGeometry::GridView::dimensionworld;

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    problem->applyInitialSolution(x);

    const Scalar Ra = getParam<Scalar>("DimensionlessNumbers.Ra");

    // Optionally start from the analytical diffusive profile at t = tStart.
    // Default: tStart = 146/Ra^2 (dimensionless critical time from linear stability).
    // Set TimeLoop.TStart = 0 in params.input to use the standard zero-concentration IC.
    const Scalar tcrit  = 146.0 / Ra;
    const Scalar tStart = getParam<Scalar>("TimeLoop.TStart", tcrit);
    const Scalar perturbAmplitude = getParam<Scalar>("Problem.PerturbAmplitude", 0.0);
    const int    perturbSeed      = getParam<int>("Problem.PerturbSeed", 1);
    if (tStart > 0.0)
    {
        const Scalar zMax = gridGeometry->bBoxMax()[dimWorld-1];

        // One N(0,1) draw per x-coordinate, memoized so every vertex in the
        // same z-column (same x, structured grid) reuses the identical draw.
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
            x[idx][Indices::concentrationIdx] = C;
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
    // Sh(t) = (1/A_top) * integral_top (dC/dn) dA  where n points into the domain
    // For pure diffusion: Sh = sqrt(Ra/(pi*t)); convection enhances this.
    const Scalar zMax = gridGeometry->bBoxMax()[dimWorld-1];

    // horizontal area of the top face (all dimensions except the vertical one)
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
                totalFlux += (grads[Indices::concentrationIdx] * scvf.unitOuterNormal()) * scvf.area();
            }
        }
        // sum over all MPI ranks
        totalFlux = leafGridView.comm().sum(totalFlux);
        return totalFlux / domainTopArea;
    };

    // --- Mass-balance flux: integrate C over domain, use dM/dt as dissolution rate ---
    auto computeTotalMass = [&]() -> Scalar {
        Scalar mass = 0.0;
        for (const auto& element : elements(leafGridView, Dune::Partitions::interior))
        {
            auto fvGeometry = localView(*gridGeometry);
            fvGeometry.bindElement(element);
            for (const auto& scv : scvs(fvGeometry))
                mass += x[scv.dofIndex()][Indices::concentrationIdx] * scv.volume();
        }
        return leafGridView.comm().sum(mass);
    };

    // --- Dominant finger count n_hat(t): energy-weighted mean wavenumber of the
    // depth-integrated vorticity field omega = -dC/dx (Riaz et al. 2006 eq. 3.10/3.11),
    // converted from radian wavenumber to a mode/finger count via /(2*pi) (valid
    // since the domain width here is 1, matching their aspect ratio A=1). ---
    constexpr Scalar pi = 3.14159265358979323846;
    const Scalar xMin = gridGeometry->bBoxMin()[0];
    const Scalar xMax = gridGeometry->bBoxMax()[0];
    const Scalar Lx   = xMax - xMin;
    int nx = 0;
    for (const auto c : getParam<std::vector<int>>("Grid.Cells0"))
        nx += c;
    const Scalar dxCol = Lx / nx;

    // Returns {n_hat, totalVorticityEnergy}. The energy is reported alongside
    // n_hat because n_hat (a ratio) is meaningless when there is no real
    // vorticity signal yet -- e.g. at PerturbAmplitude=0 the field has no
    // x-dependence in exact arithmetic, but Newton/linear-solver residual noise
    // (see Newton.MaxRelativeShift) breaks that symmetry at the ~1e-5..1e-6
    // level and can otherwise masquerade as an arbitrary "dominant" mode.
    // Compare totalVorticityEnergy against its own early-time value to judge
    // whether n_hat reflects a real, growing perturbation or just solver noise.
    auto computeDominantWavenumberAndEnergy = [&]() -> std::pair<Scalar, Scalar> {
        // depth-integrate omega into nx x-columns via a per-element midpoint/rectangle rule
        std::vector<Scalar> g(nx, 0.0);
        for (const auto& element : elements(leafGridView, Dune::Partitions::interior))
        {
            auto fvGeometry = localView(*gridGeometry);
            fvGeometry.bindElement(element);
            const auto elemSol = elementSolution(element, x, *gridGeometry);
            const auto geo = element.geometry();
            const auto center = geo.center();

            const auto grads = evalGradients(element, geo, *gridGeometry, elemSol, center);
            const Scalar omega = -grads[Indices::concentrationIdx][0]; // omega = -dC/dx

            int col = static_cast<int>(std::floor((center[0] - xMin) / dxCol));
            col = std::min(std::max(col, 0), nx - 1);
            g[col] += omega * geo.volume(); // area-integral; /dxCol below turns this into a depth-integral
        }
        for (auto& v : g)
            v = leafGridView.comm().sum(v) / dxCol;

        // energy-weighted mean wavenumber up to the Nyquist mode
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

    // --- Finger tip depth: zMax minus the shallowest z at which C >= threshold
    // for any x (i.e. how far the most advanced finger has descended). The
    // threshold is not reported by the paper -- chosen and documented here. ---
    const Scalar tipThreshold = getParam<Scalar>("Problem.TipThreshold", 0.01);
    auto computeTipDepth = [&]() -> Scalar {
        Scalar minZ = zMax;
        for (const auto& vertex : vertices(leafGridView))
        {
            const auto idx = gridGeometry->dofMapper().index(vertex);
            if (x[idx][Indices::concentrationIdx] >= tipThreshold)
                minZ = std::min(minZ, vertex.geometry().center()[dimWorld-1]);
        }
        minZ = leafGridView.comm().min(minZ);
        return zMax - minZ;
    };

    std::ofstream shFile(problem->name() + "_sherwood.csv");
    shFile << std::scientific << std::setprecision(6);
    // Sh_grad / F_grad: gradient at top interface (diffusive flux proxy)
    // Sh_mass / F_mass: dM/dt mass balance (total flux, diffusion + advection)
    shFile << "time,Sh_grad,F_grad,Sh_mass,F_mass\n";
    // -------------------------------------------------------------------------

    std::ofstream fingerFile(problem->name() + "_finger.csv");
    fingerFile << std::scientific << std::setprecision(6);
    // n_hat: dominant finger count (energy-weighted mean wavenumber / 2*pi) --
    //   only meaningful once vorticity_energy is well above its early-time
    //   (noise-floor) value, see computeDominantWavenumberAndEnergy above.
    // tip_depth: zMax - shallowest z with C >= Problem.TipThreshold (default 0.01)
    fingerFile << "time,n_hat,vorticity_energy,tip_depth\n";

    const auto tEnd       = getParam<Scalar>("TimeLoop.TEnd");
    const auto dtInit     = getParam<Scalar>("TimeLoop.DtInitial");
    const auto maxDt      = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    const auto vtkInterval = getParam<Scalar>("Output.VtkOutputInterval", 10.0);

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

    using LinearSolver = Dumux::UMFPackIstlSolver<SeqLinearSolverTraits,
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