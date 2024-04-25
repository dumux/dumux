// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPTests
 * \brief Test for the two-phase porousmedium flow model.
 */
#include <config.h>

#include <iostream>
#include <sstream>

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
#include <dumux/io/grid/gridmanager_sub.hh>
#include <dumux/io/loadsolution.hh>

#include <dumux/porousmediumflow/droplet/dropsolver.hh>
#include <dumux/porousmediumflow/droplet/newtonsolverdroplet.hh>
#include <dumux/porousmediumflow/droplet/timeloopdroplet.hh>
#include "properties.hh"

#ifndef DIFFMETHOD
#define DIFFMETHOD DiffMethod::numeric
#endif


/*!
 * \brief A method providing an () operator in order to select elements for a subgrid.
 */
template<class GlobalPosition, class Scalar>
class CircleSelector
{
public:
    CircleSelector(const GlobalPosition center, Scalar radius)
    : center_(center)
    , radius_(radius)
    {}

    //! Select all elements within a circle around a center point.
    template<class Element>
    bool operator() (const Element& element) const
    {
        const auto x = element.geometry().center()[0];
        const auto z = element.geometry().center()[2];
        // const double radius = 0.3;
        return std::hypot(x-center_[0], z-center_[2]) < radius_;
    }
private:
    const GlobalPosition center_;
    const Scalar radius_;
};

int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = Properties::TTag::TYPETAG;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    // try to create a grid (from the given grid file or the input file)
    constexpr int dim = 3;
    using HostGrid = Dune::YaspGrid<dim>;
    GridManager<HostGrid> hostGridManager;
    hostGridManager.init();
    auto& hostGrid = hostGridManager.grid();

    const auto& hostGridLeafGridView = hostGridManager.grid().leafGridView();

    // Calculate the bounding box of the host grid view.
    using GlobalPosition = Dune::FieldVector<double, dim>;
    GlobalPosition bBoxMin(std::numeric_limits<double>::max());
    GlobalPosition bBoxMax(std::numeric_limits<double>::min());
    for (const auto& vertex : vertices(hostGrid.leafGridView()))
    {
        for (int i=0; i<dim; i++)
        {
            using std::min;
            using std::max;
            bBoxMin[i] = min(bBoxMin[i], vertex.geometry().corner(0)[i]);
            bBoxMax[i] = max(bBoxMax[i], vertex.geometry().corner(0)[i]);
        }
    }
    GlobalPosition tabletCenter{0.0};
    for (int i=0; i<dim; i++)
        tabletCenter[i] = bBoxMin[i]+0.5*bBoxMax[i];

    Scalar tabletRadius = getParam<Scalar>("Tablet.Radius", 1e-3);
    CircleSelector<GlobalPosition, Scalar> elementSelector(tabletCenter, tabletRadius);
    //using SubGrid = Dune::SubGrid<3, HostGrid>;
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init(hostGrid, elementSelector);


    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    // the problem (initial and boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry, tabletCenter, tabletRadius);

    // get some time loop parameters
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");
    const auto inkjetPrintingFrequency = getParam<Scalar>("TimeLoop.InkJetFrequency");

    // check if we are about to restart a previously interrupted simulation
    Scalar restartTime = getParam<Scalar>("Restart.Time", 0);


    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());

    using DropSolver = DropletSolverTwoP<TypeTag, false>;
    auto dropSolver = std::make_shared<DropSolver>(*problem, x);
    problem->setDropSolver(dropSolver);

    problem->applyInitialSolution(x);
    auto xOld = x;

    // maybe update the interface parameters
    if (ENABLEINTERFACESOLVER)
        problem->spatialParams().updateMaterialInterfaces(x);

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);
    dropSolver->setGridVariables(gridVariables);

    // initialize the vtk output module
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;

    // use non-conforming output for the test with interface solver
    const auto ncOutput = getParam<bool>("Problem.UseNonConformingOutput", false);
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name(), "",
                                                             ncOutput ? Dune::VTK::nonconforming : Dune::VTK::conforming);
    using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.write(restartTime);

    // instantiate time loop
    auto timeLoop = std::make_shared<TimeLoopDroplet<Scalar>>(restartTime, dt, tEnd, inkjetPrintingFrequency);
    timeLoop->setMaxTimeStepSize(maxDt);
    dropSolver->setTimeLoop(timeLoop);

    // the assembler with time loop for instationary problem
    using Assembler = FVAssembler<TypeTag, DIFFMETHOD>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);

    // the linear solver
    using LinearSolver = ILURestartedGMResIstlSolver<LinearSolverTraits<GridGeometry>, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());

    // the non-linear solver
    using NewtonSolver = Dumux::NewtonSolverDroplet<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    // time loop
    timeLoop->start(); do
    {
        // solve the non-linear system with time step control
        nonLinearSolver.solve(x, *timeLoop);

        // make the new solution the old solution
        xOld = x;
        gridVariables->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        dropSolver->update();
        // write vtk output
        vtkWriter.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // Dispense a droplet at the new time
        dropSolver->dispenseDroplet();

        // set new dt as suggested by the Newton solver
        Scalar suggestedTimeStepSize = std::min(dropSolver->suggestTimeStepSize(), nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

        timeLoop->setTimeStepSize(suggestedTimeStepSize);

    } while (!timeLoop->finished());

    // output some Newton statistics
    nonLinearSolver.report();

    timeLoop->finalize(leafGridView.comm());

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
