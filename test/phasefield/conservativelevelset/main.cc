// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ConservativeLevelSetTests
 * \brief Test for the conservative level-set model discretized with PQ2
 *        (hybrid CVFE) elements.
 */
#include <config.h>

#include <cmath>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/timeloop.hh>

#include <dumux/discretization/cvfe/interpolate.hh>
#include <dumux/phasefield/common/integrate.hh>

#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/assembly/assembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dune/grid/yaspgrid.hh>
#include <dumux/io/chrono.hh>
#include <dumux/io/format.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#if DUMUX_HAVE_GRIDFORMAT
#include <dumux/io/gridwriter.hh>
#endif

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    // We initialize MPI and/or multithreading backend.
    Dumux::initialize(argc, argv);

    // parse the parameter tree (params.input)
    Parameters::init(argc, argv);

    using TypeTag = Properties::TTag::TestConservativeLevelSetPQ2;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;

    // create the grid, grid geometry, and problem
    GridManager<Grid> gridManager;
    gridManager.init();

    auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());
    auto problem = std::make_shared<Problem>(gridGeometry);

    // create the initial solution: a circular interface with a deliberately
    // SHARPER profile than the model's own equilibrium width
    // (Problem.InterfaceThickness). The reinitialization equation should then
    // smooth the profile toward its fixed-width equilibrium shape while
    // leaving the interface location (and hence the enclosed area) essentially
    // unchanged -- the classic Olsson-Kreiss reinitialization benchmark.
    SolutionVector sol(gridGeometry->numDofs());
    const auto radius = getParam<Scalar>("Problem.InitialRadius", 0.25);
    const auto interfaceWidth = getParam<Scalar>("Problem.InitialInterfaceWidth", 0.005);
    Dumux::CVFE::interpolate(*gridGeometry, sol, [&](const auto& globalPos)
    {
        using std::tanh; using std::sqrt;
        const auto dx = globalPos[0] - 0.5;
        const auto dy = globalPos[1] - 0.5;
        const auto distanceToCenter = sqrt(dx*dx + dy*dy);

        PrimaryVariables values(0.0);
        values[Indices::phaseFieldIdx] = 0.5*(1.0 - tanh((distanceToCenter - radius)/(sqrt(2.0)*interfaceWidth)));
        return values;
    });

    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(sol);

    // The phase field is a smoothed indicator of the enclosed region, so its
    // domain integral approximates the enclosed area. This equation is in the
    // same conservative (divergence) form as Cahn-Hilliard's mass balance, so
    // the quantity the CV flux scheme conserves exactly (up to Newton
    // tolerance) is likewise the control-volume (vertex) dofs' own lumped
    // sum, not the full hybrid field (see PhaseField::integrateControlVolumeField).
    auto enclosedArea = [&](const SolutionVector& s)
    {
        return PhaseField::integrateControlVolumeField(*gridGeometry, gridVariables->curGridVars(), s,
            [](const auto& vars) { return vars.phaseField(); });
    };
    const auto initialArea = enclosedArea(sol);

    // instantiate the time loop using start/end (pseudo-)time and initial time step size from params.input
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(
        Chrono::toSeconds(getParam("TimeLoop.TStart", "0")),
        Chrono::toSeconds(getParam("TimeLoop.InitialTimeStepSize")),
        Chrono::toSeconds(getParam("TimeLoop.TEnd"))
    );
    timeLoop->setMaxTimeStepSize(Chrono::toSeconds(getParam("TimeLoop.MaxTimeStepSize")));

    // assembler, linear solver, and Newton solver
    using Assembler = Dumux::Experimental::Assembler<TypeTag, DiffMethod::numeric>;
    using LinearSolver = SSORBiCGSTABIstlSolver<LinearSolverTraits<GridGeometry>, LinearAlgebraTraitsFromAssembler<Assembler>>;
    using Solver = Dumux::NewtonSolver<Assembler, LinearSolver>;

    auto oldSol = sol;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, oldSol);
    auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());
    Solver solver(assembler, linearSolver);

    // write the phase-field, at the actual (order-2) Lagrange resolution of
    // the PQ2 discretization, using GridFormat.
    // Output.Every controls how many timesteps pass between writes (1 = every step).
#if DUMUX_HAVE_GRIDFORMAT
    const int outputEvery = getParam<int>("Output.Every", 1);
    std::vector<Scalar> phaseField(sol.size());
    IO::GridWriter<typename GridGeometry::GridView, 2> vtkWriter{IO::Format::vtu, gridGeometry->gridView(), IO::order<2>};
    auto writeVtk = [&]()
    {
        for (std::size_t i = 0; i < sol.size(); ++i)
            phaseField[i] = sol[i][Indices::phaseFieldIdx];

        vtkWriter.setPointField("phi", phaseField);
        vtkWriter.write(problem->name() + "-" + Fmt::format("{:05d}", timeLoop->timeStepIndex()));
    };
    writeVtk();
#endif

    // the time loop
    timeLoop->start(); do
    {
        solver.solve(sol, *timeLoop);

        oldSol = sol;
        gridVariables->advanceTimeStep();

        timeLoop->advanceTimeStep();
        timeLoop->reportTimeStep();
        timeLoop->setTimeStepSize(solver.suggestTimeStepSize(timeLoop->timeStepSize()));

#if DUMUX_HAVE_GRIDFORMAT
        if (timeLoop->timeStepIndex() % outputEvery == 0 || timeLoop->finished())
            writeVtk();
#endif

    } while (!timeLoop->finished());
    timeLoop->finalize(gridGeometry->gridView().comm());

    const auto finalArea = enclosedArea(sol);
    const auto relativeAreaDrift = std::abs(finalArea - initialArea)/initialArea;
    std::cout << "Area conservation check: initial = " << initialArea << ", final = " << finalArea
              << ", relative drift = " << relativeAreaDrift << std::endl;

    static constexpr Scalar relativeAreaTolerance = 1e-8;
    if (relativeAreaDrift > relativeAreaTolerance)
        DUNE_THROW(Dune::Exception, "Conservative level-set enclosed area not preserved: relative drift = "
                                     << relativeAreaDrift << " > tolerance " << relativeAreaTolerance);

    return 0;
}
