// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Stationary test for the staggered grid Navier-Stokes model (Kovasznay 1948, \cite Kovasznay1948).
 */

#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/istl/io.hh>

#include <dumux/assembly/staggeredfvassembler.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/io/grid/gridmanager_sub.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/staggeredvtkoutputmodule.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>

#include <test/freeflow/navierstokes/analyticalsolutionvectors.hh>
#include <test/freeflow/navierstokes/errors.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = Properties::TTag::KovasznayTest;

    // maybe initialize MPI and/or multithreading backend
    initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // create a grid
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    Dumux::GridManager<Grid> gridManager;

#if HAVE_DUNE_SUBGRID
    const bool isStaircaseGeometry = getParam<bool>("Problem.IsStaircaseGeometry", false);

    // cut out elements within the stair-case region
    auto selector = [&](const auto& element)
    {
        if (!isStaircaseGeometry)
            return true;

        const auto globalPos = element.geometry().center();
        const auto x = globalPos[0];
        const auto y = globalPos[1];

        auto eps = 1e-12;
        static const auto lowerLeft = getParam<std::decay_t<decltype(globalPos)>>("Grid.LowerLeft");

        if (x < eps || y > (lowerLeft[1] + 1.0 - eps))
            return true;

        const bool isBelowCurve = y < (x/2.0 - std::floor(x/2.0) + lowerLeft[1] + eps);
        return !isBelowCurve;
    };

    gridManager.init(selector, "Internal");
#else
    gridManager.init();
#endif

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    // the problem (boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x;
    x[GridGeometry::cellCenterIdx()].resize(gridGeometry->numCellCenterDofs());
    x[GridGeometry::faceIdx()].resize(gridGeometry->numFaceDofs());

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // initialize the vtk output module
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    StaggeredVtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields

    NavierStokesAnalyticalSolutionVectors analyticalSolVectors(problem);
    vtkWriter.addField(analyticalSolVectors.getAnalyticalPressureSolution(), "pressureExact");
    vtkWriter.addField(analyticalSolVectors.getAnalyticalVelocitySolution(), "velocityExact");
    vtkWriter.addFaceField(analyticalSolVectors.getAnalyticalVelocitySolutionAtDofs(), "faceVelocityExact");

    vtkWriter.write(0.0);

    // the assembler with time loop for instationary problem
    using Assembler = StaggeredFVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    // the linear solver
    using LinearSolver = Dumux::UzawaBiCGSTABIstlSolver<
        SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>
    >; // if rho \neq 1, use UMFPack rather than BIGCStab
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    // linearize & solve
    Dune::Timer timer;
    nonLinearSolver.solve(x);

    // print discrete L2 and Linfity errors
    const bool printErrors = getParam<bool>("Problem.PrintErrors", false);
    const bool printConvergenceTestFile = getParam<bool>("Problem.PrintConvergenceTestFile", false);

    if (printErrors || printConvergenceTestFile)
    {
        NavierStokesErrors errors(problem, x);
        NavierStokesErrorCSVWriter(
            problem, std::to_string(x[GridGeometry::cellCenterIdx()].size())
        ).printErrors(errors);

        if (printConvergenceTestFile)
            convergenceTestAppendErrors(problem, errors);
    }

    // write vtk output
    vtkWriter.write(1.0);

    timer.stop();

    const auto& comm = Dune::MPIHelper::getCommunication();
    std::cout << "Simulation took " << timer.elapsed() << " seconds on "
              << comm.size() << " processes.\n"
              << "The cumulative CPU time was " << timer.elapsed()*comm.size() << " seconds.\n";


    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
} // end main
