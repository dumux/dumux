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
 *
 * 2D (default): parallel, h-adaptive (refine on grad(C)), optionally dynamically
 * load-balanced, on ALUGrid<2,2,simplex,conforming> -- see ../NOTES.md for the
 * development/debugging history of this machinery. The former static YaspGrid
 * boundary-layer grading is replaced by adaptive refinement (Adaptive.MaxLevel), and
 * the linear solver is the parallel ILU+GMRES this test suite standardized on (AMG
 * was investigated and rejected, see ../NOTES.md).
 *
 * 3D (SERIAL_YASP_3D compile definition, TYPETAG=BoussinesqOneSidedRB3D): unchanged
 * serial YaspGrid + UMFPack path -- the 3D conforming-bisection ALUGrid variant is
 * untested territory for the rebalance machinery and deliberately out of scope.
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

#include <dumux/discretization/evalgradients.hh>
#include <dumux/discretization/cvfe/elementsolution.hh>

#if SERIAL_YASP_3D
#include <dumux/io/grid/gridmanager_yasp.hh>
#else
#include <dumux/io/grid/gridmanager_alu.hh>

#include <dumux/adaptive/adapt.hh>
#include <dumux/adaptive/markelements.hh>
#include <dumux/adaptive/initializationindicator.hh>

#include "../common/adaptive/gridadaptindicator.hh"
#include "../common/adaptive/griddatatransfer.hh"
#include "../common/adaptive/boxrebalance.hh"
#endif

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

#if !SERIAL_YASP_3D
    // Optional vertical-strip partitioning (one x-interval per rank, full column top to
    // bottom, see loadbalancer.hh's VerticalStripDestinations): fingers grow vertically and
    // the refined boundary layer spans the full width, so strips spread both evenly over
    // the ranks and a finger stays on its rank for its whole descent. Applied here to the
    // initial partition (equal-width on the unrefined grid) and below on every rebalance
    // (weighted, preserving the strip layout).
    const bool stripPartitioning = getParam<bool>("Adaptive.VerticalStripPartitioning", false);
    if (stripPartitioning && gridManager.grid().comm().size() > 1)
    {
        BoussinesqAdaptive::VerticalStripDestinations<GetPropType<TypeTag, Properties::Grid>>
            destinations(gridManager.grid());
        gridManager.grid().repartition(destinations);
    }
#endif

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
    // Default: tStart = 146/Ra (dimensionless critical time from linear stability).
    // Set TimeLoop.TStart = 0 in params.input to use the standard zero-concentration IC.
    const Scalar tcrit  = 146.0 / Ra;
    const Scalar tStart = getParam<Scalar>("TimeLoop.TStart", tcrit);
    const Scalar perturbAmplitude = getParam<Scalar>("Problem.PerturbAmplitude", 0.0);
    const int    perturbSeed      = getParam<int>("Problem.PerturbSeed", 1);

    // One N(0,1) draw per x-coordinate, memoized so every vertex in the same z-column
    // reuses the identical draw. Reapplied in full after every initial-refinement step:
    // the rng is reseeded per application, so the noise is deterministic per run
    // configuration.
    std::size_t numClipped = 0, numVertices = 0;
    const auto setInitialConcentration = [&](SolutionVector& sol)
    {
        if (tStart <= 0.0)
            return;

        const Scalar zMax = gridGeometry->bBoxMax()[dimWorld-1];

        std::mt19937 rng(static_cast<unsigned int>(perturbSeed));
        std::normal_distribution<Scalar> noiseDist(0.0, 1.0);
        std::map<Scalar, Scalar> noiseByX;
        const auto noiseAt = [&](Scalar xCoord) -> Scalar {
            auto it = noiseByX.find(xCoord);
            if (it == noiseByX.end())
                it = noiseByX.emplace(xCoord, noiseDist(rng)).first;
            return it->second;
        };

        numClipped = 0;
        numVertices = 0;
        for (const auto& vertex : vertices(gridGeometry->gridView()))
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
            sol[idx][Indices::concentrationIdx] = C;
        }
    };
    setInitialConcentration(x);
    auto xOld = x;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

#if !SERIAL_YASP_3D
    // --- initial h-adaptive refinement (replaces the former static y-grading) ---
    // Absolute-gradient criterion (see BoussinesqGradientMagnitudeIndicator): refine where
    // the true |grad C| exceeds refineBound, coarsen where it falls below coarsenBound --
    // physical thresholds scaled by Ra, not fractions of the field's current range. Default
    // coarsenBound (1/Ra) is larger than the default refineBound (1/(5*Ra)): in the resulting
    // overlap, the coarsen check wins (see the class doc) -- only matters for the (rare, this
    // is a steep-front/near-zero bimodal field) case where a gradient magnitude falls between
    // the two.
    const Scalar refineGradBound  = getParam<Scalar>("Adaptive.RefineGradientThreshold", 1.0/(5.0*Ra));
    const Scalar coarsenGradBound = getParam<Scalar>("Adaptive.CoarsenGradientThreshold", 1.0/Ra);
    const auto maxLevel     = getParam<std::size_t>("Adaptive.MaxLevel", 0);

    BoussinesqGradientMagnitudeIndicator<TypeTag> indicator(gridGeometry, Indices::concentrationIdx);
    BoussinesqBoxGridDataTransfer<TypeTag> dataTransfer(gridGeometry, x, xOld);

    GridAdaptInitializationIndicator<TypeTag> initIndicator(problem, gridGeometry, gridVariables);
    for (std::size_t i = 0; i < maxLevel; ++i)
    {
        initIndicator.calculate(x);
        bool wasAdapted = false;
        if (markElements(gridManager.grid(), initIndicator))
            wasAdapted = adapt(gridManager.grid(), dataTransfer);

        if (wasAdapted)
        {
            setInitialConcentration(x);
            xOld = x;
            gridVariables->updateAfterGridAdaption(x);
        }
    }

    for (std::size_t i = 0; i < maxLevel; ++i)
    {
        indicator.calculate(x, refineGradBound, coarsenGradBound);
        bool wasAdapted = false;
        if (markElements(gridManager.grid(), indicator))
            wasAdapted = adapt(gridManager.grid(), dataTransfer);

        if (wasAdapted)
        {
            setInitialConcentration(x);
            xOld = x;
            gridVariables->updateAfterGridAdaption(x);
        }
        else
            break;
    }

    // The loop above's last indicator.calculate() call (at the top of its final iteration)
    // was computed *before* that iteration's adapt() -- so indicator.values() is still sized
    // for the pre-adapt mesh if the last iteration did adapt. Refresh once more so it's valid
    // for the actual final mesh before adaptIndicatorValue/adaptMark are registered with the
    // VTK writer below (a mismatch here throws Dune::RangeError at addField()).
    indicator.calculate(x, refineGradBound, coarsenGradBound);
#endif

    if (perturbAmplitude > 0.0 && tStart > 0.0 && mpiHelper.rank() == 0)
        std::cout << "Perturbed IC: amplitude = " << perturbAmplitude
                  << ", seed = " << perturbSeed
                  << ", clipped " << numClipped << "/" << numVertices
                  << " vertices to [0,1]\n";

    // --- Sherwood number: dimensionless dissolution flux at the top boundary ---
    // Sh(t) = (1/A_top) * integral_top (dC/dn) dA  where n points into the domain
    // For pure diffusion: Sh = sqrt(Ra/(pi*t)); convection enhances this.
    // All diagnostics iterate gridGeometry->gridView() (never the leafGridView captured
    // above): after adapt()/loadBalance() the initial view handle is stale on ALUGrid.
    const Scalar zMax = gridGeometry->bBoxMax()[dimWorld-1];

    // horizontal area of the top face (all dimensions except the vertical one)
    Scalar domainTopArea = 1.0;
    for (int d = 0; d < dimWorld-1; ++d)
        domainTopArea *= gridGeometry->bBoxMax()[d] - gridGeometry->bBoxMin()[d];

    auto computeSherwood = [&]() -> Scalar {
        Scalar totalFlux = 0.0;
        for (const auto& element : elements(gridGeometry->gridView(), Dune::Partitions::interior))
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
        totalFlux = gridGeometry->gridView().comm().sum(totalFlux);
        return totalFlux / domainTopArea;
    };

    // --- Mass-balance flux: integrate C over domain, use dM/dt as dissolution rate ---
    auto computeTotalMass = [&]() -> Scalar {
        Scalar mass = 0.0;
        for (const auto& element : elements(gridGeometry->gridView(), Dune::Partitions::interior))
        {
            auto fvGeometry = localView(*gridGeometry);
            fvGeometry.bindElement(element);
            for (const auto& scv : scvs(fvGeometry))
                mass += x[scv.dofIndex()][Indices::concentrationIdx] * scv.volume();
        }
        return gridGeometry->gridView().comm().sum(mass);
    };

    // --- Dominant finger count n_hat(t): energy-weighted mean wavenumber of the
    // depth-integrated vorticity field omega = -dC/dx (Riaz et al. 2006 eq. 3.10/3.11),
    // converted from radian wavenumber to a mode/finger count via /(2*pi) (valid
    // since the domain width here is 1, matching their aspect ratio A=1). ---
    constexpr Scalar pi = 3.14159265358979323846;
    const Scalar xMin = gridGeometry->bBoxMin()[0];
    const Scalar xMax = gridGeometry->bBoxMax()[0];
    const Scalar Lx   = xMax - xMin;
    // Fixed x-column binning resolution for the depth integration: sum of the YaspGrid
    // Cells0 blocks if given (legacy key, still used by the pressure test sharing this
    // params file), else the first entry of the ALUGrid Cells key. The binning is a
    // diagnostic resolution only and stays fixed while the grid adapts underneath.
    int nx = 0;
    if (hasParam("Grid.Cells0"))
        for (const auto c : getParam<std::vector<int>>("Grid.Cells0"))
            nx += c;
    else
        nx = getParam<std::vector<int>>("Grid.Cells")[0];
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
        for (const auto& element : elements(gridGeometry->gridView(), Dune::Partitions::interior))
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
            v = gridGeometry->gridView().comm().sum(v) / dxCol;

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
        for (const auto& vertex : vertices(gridGeometry->gridView()))
        {
            const auto idx = gridGeometry->dofMapper().index(vertex);
            if (x[idx][Indices::concentrationIdx] >= tipThreshold)
                minZ = std::min(minZ, vertex.geometry().center()[dimWorld-1]);
        }
        minZ = gridGeometry->gridView().comm().min(minZ);
        return zMax - minZ;
    };

    // CSV output on rank 0 only (every rank still participates in the comm()
    // reductions inside the compute functions above).
    const bool isIORank = (mpiHelper.rank() == 0);
    std::ofstream shFile, fingerFile;
    if (isIORank)
    {
        shFile.open(problem->name() + "_sherwood.csv");
        shFile << std::scientific << std::setprecision(6);
        // Sh_grad / F_grad: gradient at top interface (diffusive flux proxy)
        // Sh_mass / F_mass: dM/dt mass balance (total flux, diffusion + advection)
        shFile << "time,Sh_grad,F_grad,Sh_mass,F_mass\n";

        fingerFile.open(problem->name() + "_finger.csv");
        fingerFile << std::scientific << std::setprecision(6);
        // n_hat: dominant finger count (energy-weighted mean wavenumber / 2*pi) --
        //   only meaningful once vorticity_energy is well above its early-time
        //   (noise-floor) value, see computeDominantWavenumberAndEnergy above.
        // tip_depth: zMax - shallowest z with C >= Problem.TipThreshold (default 0.01)
        fingerFile << "time,n_hat,vorticity_energy,tip_depth\n";
    }

    const auto tEnd       = getParam<Scalar>("TimeLoop.TEnd");
    const auto dtInit     = getParam<Scalar>("TimeLoop.DtInitial");
    const auto maxDt      = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    const auto vtkInterval = getParam<Scalar>("Output.VtkOutputInterval", 10.0);

    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    IOFields::initOutputModule(vtkWriter);

#if !SERIAL_YASP_3D
    // --- adaptivity diagnostics ---
    // "adaptIndicatorValue": indicator.values(), the true |grad C| at each element's center
    // that operator() compares against refineBound()/coarsenBound() (see
    // gridadaptindicator.hh's BoussinesqGradientMagnitudeIndicator).
    // Passed by reference once: the same std::vector object is what indicator.calculate()
    // updates in place every timestep, so this field tracks it automatically -- *except*
    // right after an adapt()/rebalance changes element indices, which is why the calculate()
    // calls below are duplicated after those events (indicator.values() must be resized to
    // match the *current* mesh before the next write(), not just the next timestep).
    // "adaptMark": operator()'s decision for each element (1 refine / -1 coarsen / 0 keep),
    // i.e. which side of the two bounds each element actually fell on. Recomputed into its
    // own persistent vector right before every write() (cheap, and operator() itself doesn't
    // change any state, so recomputing it doesn't need the calculate() re-run above).
    std::vector<Scalar> adaptMarkField;
    const auto updateAdaptMarkField = [&]() {
        const auto& gv = gridGeometry->gridView();
        adaptMarkField.assign(gv.size(0), 0.0);
        for (const auto& element : elements(gv))
            adaptMarkField[gridGeometry->elementMapper().index(element)] = static_cast<Scalar>(indicator(element));
    };
    updateAdaptMarkField();
    vtkWriter.addField(indicator.values(), "adaptIndicatorValue", Vtk::FieldType::element);
    vtkWriter.addField(adaptMarkField, "adaptMark", Vtk::FieldType::element);
#endif

    vtkWriter.write(0.0);

    const auto writeSherwood = [&](Scalar t, Scalar massFluxPerWidth = 0.0) {
        const Scalar Sh_grad = computeSherwood();
        const Scalar Sh_mass = massFluxPerWidth * Ra;
        if (isIORank)
        {
            shFile << t << "," << Sh_grad << "," << Sh_grad/Ra
                        << "," << Sh_mass << "," << massFluxPerWidth << "\n";
            shFile.flush();
        }
    };

    const auto writeFinger = [&](Scalar t) {
        const auto [nHat, vorticityEnergy] = computeDominantWavenumberAndEnergy();
        const Scalar tipDepth = computeTipDepth();
        if (isIORank)
        {
            fingerFile << t << "," << nHat << "," << vorticityEnergy << "," << tipDepth << "\n";
            fingerFile.flush();
        }
    };

    Scalar massOld = computeTotalMass();
    writeSherwood(tStart);
    writeFinger(tStart);

    Scalar nextVtkTime = tStart + vtkInterval;
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(tStart, dtInit, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);

#if SERIAL_YASP_3D
    using LinearSolver = Dumux::UMFPackIstlSolver<SeqLinearSolverTraits,
                         LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();
#else
    // Parallel ILU+GMRES: the solver the adaptive tests standardized on (AMG rejected,
    // see ../NOTES.md). Constructed from gridGeometry->gridView(), which is
    // fresh after the initial refinement above.
    using LinearSolver = ILURestartedGMResIstlSolver<LinearSolverTraits<GridGeometry>,
                                                     LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());
#endif

    NewtonSolver<Assembler, LinearSolver> nonLinearSolver(assembler, linearSolver);

#if !SERIAL_YASP_3D
    // --- dynamic load balancing (opt-in), see common/adaptive/{loadbalancer,boxrebalance}.hh ---
    const bool enableRebalancing = getParam<bool>("Adaptive.EnableDynamicLoadBalancing", false);
    const auto rebalanceEvery = getParam<std::size_t>("Adaptive.RebalanceEvery", 5);
    const Scalar imbalanceTolerance = getParam<Scalar>("Adaptive.ImbalanceTolerance", 0.15);
    std::size_t stepsSinceRebalance = 0;

    const auto localInteriorElements = [&]() -> Scalar {
        Scalar n = 0;
        for (const auto& element : elements(gridGeometry->gridView(), Dune::Partitions::interior))
        { (void)element; n += 1; }
        return n;
    };

    const auto currentImbalance = [&]() -> Scalar {
        const Scalar local = localInteriorElements();
        const auto& comm = gridGeometry->gridView().comm();
        const Scalar maxLocal = comm.max(local);
        const Scalar sumLocal = comm.sum(local);
        const Scalar meanLocal = sumLocal / static_cast<Scalar>(comm.size());
        return meanLocal > 0.0 ? (maxLocal / meanLocal - 1.0) : 0.0;
    };
#endif

    std::size_t stepIdx = 0;
    timeLoop->start(); do
    {
#if !SERIAL_YASP_3D
        // adapt only between timesteps: at this point x == xOld (both hold the solution
        // at the end of the previous accepted timestep), so interpolating both through
        // the data transfer is consistent -- no separate "old" state is lost.
        if (stepIdx > 0)
        {
            indicator.calculate(x, refineGradBound, coarsenGradBound);

            bool wasAdapted = false;
            if (markElements(gridManager.grid(), indicator))
                wasAdapted = adapt(gridManager.grid(), dataTransfer);

            if (wasAdapted)
            {
                gridVariables->updateAfterGridAdaption(x);
                assembler->updateAfterGridAdaption();
                linearSolver->updateAfterGridAdaption(gridGeometry->gridView(), gridGeometry->dofMapper());

                // indicator.values() was sized/indexed for the pre-adapt mesh; refresh it
                // against the new one so it stays valid for adaptMarkField/VTK output below
                // (element indices shift on adapt, not just their count).
                indicator.calculate(x, refineGradBound, coarsenGradBound);
            }

            ++stepsSinceRebalance;
            if (enableRebalancing && gridGeometry->gridView().comm().size() > 1
                && (stepsSinceRebalance >= rebalanceEvery || currentImbalance() > imbalanceTolerance))
            {
                if (BoussinesqAdaptive::rebalanceBox<TypeTag>(gridGeometry, gridManager.grid(), x, xOld, stripPartitioning))
                {
                    gridVariables->updateAfterGridAdaption(x);
                    assembler->updateAfterGridAdaption();
                    linearSolver->updateAfterGridAdaption(gridGeometry->gridView(), gridGeometry->dofMapper());

                    // same reason as above: rebalancing also moves elements between (and
                    // renumbers them within) ranks
                    indicator.calculate(x, refineGradBound, coarsenGradBound);
                }
                stepsSinceRebalance = 0;
            }

            // interpolation during adapt (and, to a lesser extent, migration during
            // rebalance) changes the discrete mass integral -- rebaseline so the dM/dt
            // Sherwood diagnostic measures physics, not the grid change
            if (wasAdapted)
                massOld = computeTotalMass();
        }
#endif

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
#if !SERIAL_YASP_3D
            updateAdaptMarkField();
#endif
            vtkWriter.write(timeLoop->time());
            nextVtkTime += vtkInterval;
        }

        writeSherwood(timeLoop->time(), massFluxPerWidth);
        writeFinger(timeLoop->time());

        timeLoop->reportTimeStep();
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));
        ++stepIdx;

    } while (!timeLoop->finished());

    timeLoop->finalize(gridGeometry->gridView().comm());

    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/false);

    return 0;
}
