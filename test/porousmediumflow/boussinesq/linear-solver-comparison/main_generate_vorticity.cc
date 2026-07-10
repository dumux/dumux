// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Matrix-generation stage for the vector-potential (vorticity)
 *        Boussinesq linear-solver benchmark.
 *
 * Runs the same transient one-sided Rayleigh-Bénard setup as
 * white-noise-perturbations/main_vorticity.cc (perturbed erfc initial
 * condition, production UMFPack solver) from t = tStart to TimeLoop.TEnd, but
 * instead of writing VTK/diagnostics output it exports the assembled Jacobian
 * and residual at each time in MatrixExport.Times (default 1, 2, 3, 4) to
 * MatrixMarket files (see matrixio.hh for the naming convention).
 *
 * These files are the input to main_solve_vorticity.cc, which benchmarks
 * several ISTL linear solvers against them. Splitting the two stages means
 * the (expensive, one-off) transient run and matrix assembly only has to
 * happen once, while the solver comparison can be re-run cheaply against the
 * same matrices.
 */
#include <config.h>

#include <cmath>
#include <iostream>
#include <map>
#include <random>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/timeloop.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/fvassembler.hh>

#include <dumux/io/grid/gridmanager_yasp.hh>

#include "properties_vorticity.hh"
#include "matrixio.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTag = Properties::TTag::BoussinesqOneSidedRB;

    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    Parameters::init(argc, argv, [](auto&){}, "params.input");
    GetPropType<TypeTag, Properties::FluidSystem>::init();

    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    const auto& leafGridView = gridManager.grid().leafGridView();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    using Scalar  = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    static constexpr int dimWorld = GridGeometry::GridView::dimensionworld;

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    problem->applyInitialSolution(x);

    const Scalar Ra = getParam<Scalar>("DimensionlessNumbers.Ra");

    // Perturbed analytical diffusive profile at t = tStart -- identical recipe
    // to white-noise-perturbations/main_vorticity.cc, so the assembled
    // Jacobians develop real x-dependence (upwind-coupled transport) instead
    // of a trivial 1D erfc profile.
    const Scalar tcrit  = 146.0 / Ra;
    const Scalar tStart = getParam<Scalar>("TimeLoop.TStart", tcrit);
    const Scalar perturbAmplitude = getParam<Scalar>("Problem.PerturbAmplitude", 0.0);
    const int    perturbSeed      = getParam<int>("Problem.PerturbSeed", 1);
    {
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
                C = std::min(std::max(C, 0.0), 1.0);
            }
            x[idx][Indices::concentrationIdx] = C;
        }
    }

    auto xOld = x;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    const auto tEnd        = getParam<Scalar>("TimeLoop.TEnd");
    const auto dtInit      = getParam<Scalar>("TimeLoop.DtInitial");
    const auto maxDt       = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    const auto exportTimes = getParam<std::vector<Scalar>>("MatrixExport.Times");

    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(tStart, dtInit, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);
    timeLoop->setCheckPoint(exportTimes.begin(), exportTimes.end());

    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);

    // Production solver for the transient run itself -- matches
    // white-noise-perturbations/main_vorticity.cc. This is *not* one of the
    // solvers benchmarked by main_solve_vorticity.cc; it just needs to be
    // robust enough to reliably advance the simulation to the checkpoints.
    using LinearSolver = Dumux::UMFPackIstlSolver<SeqLinearSolverTraits,
                         LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();
    NewtonSolver<Assembler, LinearSolver> nonLinearSolver(assembler, linearSolver);

    if (mpiHelper.rank() == 0)
    {
        std::cout << "Generating vorticity-formulation matrices at t = ";
        for (const auto t : exportTimes)
            std::cout << t << " ";
        std::cout << "(TEnd = " << tEnd << ")\n";
    }

    timeLoop->start(); do
    {
        nonLinearSolver.solve(x, *timeLoop);

        xOld = x;
        gridVariables->advanceTimeStep();
        timeLoop->advanceTimeStep();

        if (timeLoop->isCheckPoint())
        {
            assembler->assembleJacobianAndResidual(x);
            exportSystem(problem->name(), timeLoop->time(), assembler->jacobian(), assembler->residual());
            if (mpiHelper.rank() == 0)
                std::cout << "Exported matrix/rhs at t = " << timeLoop->time() << "\n";
        }

        timeLoop->reportTimeStep();
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/false);

    return 0;
}
