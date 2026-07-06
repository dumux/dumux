// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup HeleShawTests
 * \brief Saffman-Taylor fingering test for the hybrid CVFE (PQ2) Hele-Shaw
 *        Darcy-Cahn-Hilliard model. No adaptivity in this first pass (fixed
 *        grid); see the classic Box variant at
 *        `test/phasefield/heleshaw/2p/` for the AMR-enabled version.
 */
#include <config.h>

#include <vector>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/timeloop.hh>

#include <dumux/discretization/cvfe/interpolate.hh>

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

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv);

    using TypeTag = Properties::TTag::TestHeleShawTwoPCVFE;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    // create the grid, grid geometry, and problem
    GridManager<Grid> gridManager;
    gridManager.init();

    auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());
    auto problem = std::make_shared<Problem>(gridGeometry);

    // create the initial solution from the problem's initialAtPos
    SolutionVector sol(gridGeometry->numDofs());
    Dumux::CVFE::interpolate(*gridGeometry, sol, [&](const auto& globalPos)
    { return problem->initialAtPos(globalPos); });

    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(sol);

    // instantiate the time loop using start/end time and initial time step size from params.input
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

    // write the pressure, phase field, and chemical potential fields, at the
    // actual (order-2) Lagrange resolution of the PQ2 discretization.
#if DUMUX_HAVE_GRIDFORMAT
    const int outputEvery = getParam<int>("Output.Every", 1);
    std::vector<Scalar> pressure(sol.size());
    std::vector<Scalar> phaseField(sol.size());
    std::vector<Scalar> chemicalPotential(sol.size());
    IO::GridWriter<typename GridGeometry::GridView, 2> vtkWriter{IO::Format::vtu, gridGeometry->gridView(), IO::order<2>};
    auto writeVtk = [&]()
    {
        for (std::size_t i = 0; i < sol.size(); ++i)
        {
            pressure[i] = sol[i][Indices::pressureIdx];
            phaseField[i] = sol[i][Indices::phaseFieldIdx];
            chemicalPotential[i] = sol[i][Indices::chemPotIdx];
        }
        vtkWriter.setPointField("p", pressure);
        vtkWriter.setPointField("phi", phaseField);
        vtkWriter.setPointField("mu", chemicalPotential);
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

    return 0;
}
