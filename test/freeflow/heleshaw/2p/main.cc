// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Saffman-Taylor fingering test for the Hele-Shaw Darcy-Cahn-Hilliard 2p model.
 *        Adaptive mesh refinement is applied to resolve the diffuse interface.
 */

#include <config.h>

#include <dune/alugrid/grid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/timeloop.hh>

#include <dumux/discretization/box.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_alu.hh>

#include <dumux/adaptive/adapt.hh>
#include <dumux/adaptive/markelements.hh>

#include <dumux/freeflow/heleshaw/2p/model.hh>
#include <dumux/freeflow/heleshaw/2p/adapt.hh>
#include "problem.hh"

// ============================================================
// TypeTag for this test
// ============================================================
namespace Dumux::Properties::TTag {

struct HeleShawTest
{
    using InheritsFrom = std::tuple<HeleShawTwoP, BoxModel>;

    using Grid = Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>;

    template<class TypeTag>
    using Problem = HeleShawTwoPTestProblem<TypeTag>;

    using EnableGridGeometryCache        = std::true_type;
    using EnableGridVolumeVariablesCache = std::true_type;
    using EnableGridFluxVariablesCache   = std::true_type;
};

} // end namespace Dumux::Properties::TTag

// ============================================================
// main
// ============================================================
int main(int argc, char** argv)
{
    using namespace Dumux;

    initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    Parameters::init(argc, argv);

    using TypeTag = Properties::TTag::HeleShawTest;

    using Scalar         = GetPropType<TypeTag, Properties::Scalar>;
    using Grid           = GetPropType<TypeTag, Properties::Grid>;
    using GridGeometry   = GetPropType<TypeTag, Properties::GridGeometry>;
    using Problem        = GetPropType<TypeTag, Properties::Problem>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables  = GetPropType<TypeTag, Properties::GridVariables>;
    using IOFields       = GetPropType<TypeTag, Properties::IOFields>;

    // ---- grid ----
    GridManager<Grid> gridManager;
    gridManager.init();

    auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());
    auto problem      = std::make_shared<Problem>(gridGeometry);

    // ---- initial solution ----
    SolutionVector sol(gridGeometry->numDofs());
    problem->applyInitialSolution(sol);
    auto oldSol = sol;

    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(sol);

    // ---- AMR setup ----
    const Scalar refineTol  = getParam<Scalar>("Adaptive.RefineTolerance");
    const Scalar coarsenTol = getParam<Scalar>("Adaptive.CoarsenTolerance");

    HeleShawTwoPGridAdaptIndicator<TypeTag> indicator(gridGeometry);
    HeleShawTwoPGridDataTransfer<TypeTag>   dataTransfer(problem, gridGeometry, gridVariables, sol);

    // Initial AMR: refine around the starting interface
    const auto initMaxLevel = getParam<std::size_t>("Adaptive.InitMaxLevel", 0);
    for (std::size_t i = 0; i < initMaxLevel; ++i)
    {
        indicator.calculate(sol, refineTol, coarsenTol);
        if (markElements(gridManager.grid(), indicator))
        {
            adapt(gridManager.grid(), dataTransfer);
            oldSol = sol;
            gridVariables->updateAfterGridAdaption(sol);
        }
    }

    // ---- VTK output ----
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, sol, problem->name());
    IOFields::initOutputModule(vtkWriter);
    vtkWriter.write(0.0);

    // ---- time loop ----
    const auto tEnd  = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt          = getParam<Scalar>("TimeLoop.DtInitial");

    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // ---- assembler + solver (rebuilt after grid adaptation) ----
    using Assembler    = FVAssembler<TypeTag, DiffMethod::numeric>;
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits,
                                            LinearAlgebraTraitsFromAssembler<Assembler>>;
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;

    auto linearSolver = std::make_shared<LinearSolver>();

    // Lambda that (re)builds assembler + Newton solver with the current old solution
    auto makeNewton = [&]() {
        auto assembler = std::make_shared<Assembler>(
            problem, gridGeometry, gridVariables, timeLoop, oldSol);
        return NewtonSolver(assembler, linearSolver);
    };

    auto newton = makeNewton();

    // ---- time loop ----
    timeLoop->start(); do
    {
        newton.solve(sol, *timeLoop);

        oldSol = sol;
        gridVariables->advanceTimeStep();
        timeLoop->advanceTimeStep();
        vtkWriter.write(timeLoop->time());
        timeLoop->reportTimeStep();
        timeLoop->setTimeStepSize(newton.suggestTimeStepSize(timeLoop->timeStepSize()));

        // AMR after each time step
        indicator.calculate(sol, refineTol, coarsenTol);
        if (markElements(gridManager.grid(), indicator))
        {
            adapt(gridManager.grid(), dataTransfer);
            oldSol = sol;
            gridVariables->updateAfterGridAdaption(sol);
            newton = makeNewton();  // rebuild for new DOF structure
        }

    } while (!timeLoop->finished());

    timeLoop->finalize(gridGeometry->gridView().comm());

    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
}
