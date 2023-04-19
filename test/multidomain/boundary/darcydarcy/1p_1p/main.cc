// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoundaryTests
 * \brief Test for the equal dimension boundary coupling model.
 */

#include <config.h>

#ifndef DOMAINSPLIT
#define DOMAINSPLIT 0
#endif

#include <iostream>
#include <fstream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/io/vtkoutputmodule.hh>

#include <dumux/io/grid/gridmanager_yasp.hh>

#if DOMAINSPLIT==1
#include <dumux/io/grid/gridmanager_sub.hh>
#endif

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/fvgridgeometry.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // Define the sub problem type tags
    using TypeTag = Properties::TTag::OnePSub;
    using SubTypeTag0 = Properties::TTag::OnePSub0;
    using SubTypeTag1 = Properties::TTag::OnePSub1;

    // create the full grid that we are gonna split for the output
    using HostGrid = typename GetProp<TypeTag, Properties::Grid>::HostGrid;
    GridManager<HostGrid> gridManager;
    gridManager.init();

    // get the lens
    using GlobalPosition = typename HostGrid::template Codim<0>::Geometry::GlobalCoordinate;
    const auto lensLowerLeft = getParam<GlobalPosition>("SpatialParams.LensLowerLeft");
    const auto lensUpperRight = getParam<GlobalPosition>("SpatialParams.LensUpperRight");

    // grid managers for the subdomains
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager0, gridManager1;

///////////////////////////////////////////////////////////////////////////////////////////
// split the domains by creating two separate grids for lens and the rest (using sub-grid)
///////////////////////////////////////////////////////////////////////////////////////////
#if DOMAINSPLIT==1

    // create subgrids coinciding with the lens
    auto elementSelector0 = [&lensLowerLeft, &lensUpperRight](const auto& element)
    { return LensSpatialParams::pointInLens(element.geometry().center(), lensLowerLeft, lensUpperRight); };
    auto elementSelector1 = [&lensLowerLeft, &lensUpperRight](const auto& element)
    { return !LensSpatialParams::pointInLens(element.geometry().center(), lensLowerLeft, lensUpperRight); };

    gridManager0.init(gridManager.grid(), elementSelector0, "1");
    gridManager1.init(gridManager.grid(), elementSelector1, "2");

//////////////////////////////////////////////////////////////////////////////////
// split the domains by creating two separate grids for the lower and upper half
//////////////////////////////////////////////////////////////////////////////////
#elif DOMAINSPLIT==0

    // create an upper half and lower half grid
    gridManager0.init("1");
    gridManager1.init("2");

#endif

    // we compute on the leaf grid views
    const auto& gridView0 = gridManager0.grid().leafGridView();
    const auto& gridView1 = gridManager1.grid().leafGridView();

    ////////////////////////////////////////////////
    // run the multidomain simulation on two grids
    ////////////////////////////////////////////////

    // the mixed dimension type traits
    using Traits = MultiDomainTraits<SubTypeTag0, SubTypeTag1>;
    constexpr auto domain0Idx = Traits::template SubDomain<0>::Index();
    constexpr auto domain1Idx = Traits::template SubDomain<1>::Index();

    // create the finite volume grid geometries
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto fvGridGeometry0 = std::make_shared<GridGeometry>(gridView0);
    auto fvGridGeometry1 = std::make_shared<GridGeometry>(gridView1);

    // the usage of multidomain fv grid geometry is unnecessarily complicated
    // but we want to test its constructor and a coupled of interfaces
    MultiDomainFVGridGeometry<Traits> mdGridGeometry(std::make_tuple(
        fvGridGeometry0, fvGridGeometry1
    ));

    // the update here is unnecessary because update is already done in the constructor
    // we do this just for testing purposes here! This means the grid geometry is updated twice.
    mdGridGeometry.update(gridView0, gridView1);

    // this is an identity operation, just to test the interface
    fvGridGeometry0 = mdGridGeometry.template get<domain0Idx>();
    fvGridGeometry1 = mdGridGeometry.template get<domain1Idx>();

    // the coupling manager
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    auto couplingManager = std::make_shared<CouplingManager>();

    // the problem (initial and boundary conditions)
    using Problem0 = GetPropType<SubTypeTag0, Properties::Problem>;
    auto problem0 = std::make_shared<Problem0>(fvGridGeometry0, couplingManager, "1");
    problem0->spatialParams().setLens(lensLowerLeft, lensUpperRight);
    using Problem1 = GetPropType<SubTypeTag1, Properties::Problem>;
    auto problem1 = std::make_shared<Problem1>(fvGridGeometry1, couplingManager, "2");
    problem1->spatialParams().setLens(lensLowerLeft, lensUpperRight);

    // the solution vector
    Traits::SolutionVector sol;
    sol[domain0Idx].resize(fvGridGeometry0->numDofs());
    sol[domain1Idx].resize(fvGridGeometry1->numDofs());
    problem0->applyInitialSolution(sol[domain0Idx]);
    problem1->applyInitialSolution(sol[domain1Idx]);
    auto oldSol = sol;

    // compute the coupling map
    couplingManager->init(problem0, problem1, sol);

    // the grid variables
    using GridVariables0 = GetPropType<SubTypeTag0, Properties::GridVariables>;
    auto gridVariables0 = std::make_shared<GridVariables0>(problem0, fvGridGeometry0);
    using GridVariables1 = GetPropType<SubTypeTag1, Properties::GridVariables>;
    auto gridVariables1 = std::make_shared<GridVariables1>(problem1, fvGridGeometry1);
    gridVariables0->init(sol[domain0Idx]);
    gridVariables1->init(sol[domain1Idx]);

    // get some time loop parameters
    using Scalar = Traits::Scalar;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // initialize the vtk output module
    using SolutionVector0 = std::decay_t<decltype(sol[domain0Idx])>;
    VtkOutputModule<GridVariables0, SolutionVector0> vtkWriter0(*gridVariables0, sol[domain0Idx], problem0->name());
    GetPropType<SubTypeTag0, Properties::IOFields>::initOutputModule(vtkWriter0);
    vtkWriter0.write(0.0);

    using SolutionVector1 = std::decay_t<decltype(sol[domain1Idx])>;
    VtkOutputModule<GridVariables1, SolutionVector1> vtkWriter1(*gridVariables1, sol[domain1Idx], problem1->name());
    GetPropType<SubTypeTag1, Properties::IOFields>::initOutputModule(vtkWriter1);
    vtkWriter1.write(0.0);

    // instantiate time loop
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // the assembler with time loop for transient problem
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(problem0, problem1),
                                                 std::make_tuple(fvGridGeometry0, fvGridGeometry1),
                                                 std::make_tuple(gridVariables0, gridVariables1),
                                                 couplingManager, timeLoop, oldSol);

    // the linear solver
    using LinearSolver = AMGBiCGSTABIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    // time loop
    timeLoop->start();
    while (!timeLoop->finished())
    {
        // solve the non-linear system with time step control
        nonLinearSolver.solve(sol, *timeLoop);

        // make the new solution the old solution
        oldSol = sol;
        gridVariables0->advanceTimeStep();
        gridVariables1->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // write vtk output
        vtkWriter0.write(timeLoop->time());
        vtkWriter1.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by newton controller
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));
    }

    timeLoop->finalize(mpiHelper.getCommunication());

    //////////////////////////////////////////////////////////////////////////
    // write out a combined vtu for vtu comparison in the testing framework
    //////////////////////////////////////////////////////////////////////////

    const auto& gridView = gridManager.grid().leafGridView();
    CCTpfaFVGridGeometry<typename HostGrid::LeafGridView> gridGeometry(gridView);
    const auto& bBoxTree = gridGeometry.boundingBoxTree();
    // copy data from the subdomains to full domain data vectors
    std::vector<int> processRank(gridView.size(0), 0); // sequential simulation
    std::vector<int> pressure(gridView.size(0));

    for (const auto& element : elements(gridView0))
    {
        const auto eIdx = fvGridGeometry0->elementMapper().index(element);
        const auto eIdxHost = intersectingEntities(element.geometry().center(), bBoxTree)[0];
        pressure[eIdxHost] = sol[domain0Idx][eIdx][0];
    }

    for (const auto& element : elements(gridView1))
    {
        const auto eIdx = fvGridGeometry1->elementMapper().index(element);
        const auto eIdxHost = intersectingEntities(element.geometry().center(), bBoxTree)[0];
        pressure[eIdxHost] = sol[domain1Idx][eIdx][0];
    }

    Dune::VTKWriter<typename HostGrid::LeafGridView> vtkWriter(gridView);
    vtkWriter.addCellData(processRank, "process rank");
    vtkWriter.addCellData(pressure, "pressure");
    const auto filename = getParam<std::string>("Vtk.OutputName") + "_combined";
    vtkWriter.write(filename);

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
}
