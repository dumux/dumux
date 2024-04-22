// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 *
 * \brief Test for the pore network model
*/
#include <config.h>
#include <ctime>
#include <iostream>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/common/properties/model.hh>
#include <dumux/common/properties/grid.hh>
#include <dumux/discretization/porenetwork/gridgeometry.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/thresholdcapillarypressures.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/pore/2p/localrulesforplatonicbody.hh>
#include <dumux/io/grid/porenetwork/gridmanager.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/transmissibility1p.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/transmissibility2p.hh>
#include <dumux/porenetwork/common/throatproperties.hh>

#include <dumux/porenetwork/2p/static/staticdrainge.hh>
#include <dumux/io/gnuplotinterface.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/pdesolver.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/assembly/fvassembler.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/porenetwork/common/pnmvtkoutputmodule.hh>

#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/grid/porenetwork/gridmanager.hh> // for pore-network grid
#include <dumux/porenetwork/common/boundaryflux.hh>

#include <dumux/porenetwork/2p/newtonsolver.hh>
#include <dumux/porenetwork/common/outletpcgradient.hh>

#include "upscalinghelper.hh"
#include "properties.hh"

// [[/codeblock]]
// [[/details]]
//
// ### The driver function
//
// The template argument `TypeTag` determines if we run the example assuming
// a creeping flow regime or not. Which regime is selected is set with the parameter
// `Problem.AssumeCreepingFlow` in the input file.
namespace Dumux {

template<class GridManager, class StaticSw>
void runSinglePhaseUpscaling(GridManager& gridManager, const StaticSw sw, int step)
{
    const auto& swStatic = sw;
    // ### Create the grid and the grid geometry
    // [[codeblock]]
    // The grid manager can be used to create a grid from the input file
    // using GridManager = PoreNetwork::GridManager<3>;
    // GridManager gridManager;
    // gridManager.init();
    using TypeTag = Properties::TTag::PNMTWOPDynamic;
    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();
    auto gridData = gridManager.getGridData();

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView, *gridData);

    // the spatial parameters
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;
    auto spatialParams = std::make_shared<SpatialParams>(gridGeometry);

    // the problem (boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry, spatialParams);
    problem->setSw(swStatic);
    problem->setName(problem->name() + std::to_string(step));

    // the solution vector
    using GridView = typename GridGeometry::GridView;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(leafGridView.size(GridView::dimension));
    problem->applyInitialSolution(x);
    auto xOld = x;

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // get some time loop parameters
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // check if we are about to restart a previously interrupted simulation
    Scalar restartTime = 0;
    if (Parameters::getTree().hasKey("Restart") || Parameters::getTree().hasKey("TimeLoop.Restart"))
        restartTime = getParam<Scalar>("TimeLoop.Restart");

    // initialize the vtk output module
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    PoreNetwork::VtkOutputModule<GridVariables, GetPropType<TypeTag, Properties::FluxVariables>, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    IOFields::initOutputModule(vtkWriter); //! Add model specific output fields

    vtkWriter.write(0.0);

    const auto outletCapPressureGradient = std::make_shared<Dumux::PoreNetwork::OutletCapPressureGradient<GridVariables, SolutionVector>>(*gridVariables, x);
    problem->outletCapPressureGradient(outletCapPressureGradient);

    // instantiate time loop
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(restartTime, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // the assembler with time loop for instationary problem
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);

    // the linear solver
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = PoreNetwork::TwoPNewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    // time loop
    timeLoop->start(); do
    {
        // try solving the non-linear system
        nonLinearSolver.solve(x, *timeLoop);

        // make the new solution the old solution
        xOld = x;
        gridVariables->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // write vtk output
        if(problem->shouldWriteOutput(timeLoop->timeStepIndex(), *gridVariables))
            vtkWriter.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by newton solver
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());

    // // ### Prepare the upscaling procedure.
    // // Set up a helper class to determine the total mass flux leaving the network
    // const auto boundaryFlux = PoreNetwork::BoundaryFlux(*gridVariables, assembler->localResidual(), x);

    // // // ### The actual upscaling procedure
    // // // #### Instantiate the upscaling helper
    // // // [[codeblock]]
    // using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    // using UpscalingHelper = UpscalingHelper<Scalar>;
    // UpscalingHelper upscalingHelper;
    // // // [[/codeblock]]

    // // // Set the side lengths used for applying the pressure gradient and calculating the REV outflow area.
    // // // One can either specify these values manually (usually more accurate) or let the UpscalingHelper struct
    // // // determine it automatically based on the network's bounding box.
    // // // [[codeblock]]
    // // const auto sideLengths = [&]()
    // // {
    // //     using GlobalPosition = typename GridGeometry::GlobalCoordinate;
    // //     if (hasParam("Problem.SideLength"))
    // //         return getParam<GlobalPosition>("Problem.SideLength");
    // //     else
    // //         return upscalingHelper.getSideLengths(*gridGeometry);
    // // }();

    // // pass the side lengths to the problem
    // // problem->setSideLengths(sideLengths);
    // // [[/codeblock]]

    // // Get the maximum and minimum pressure gradient and the population of sample points specified in the input file
    // // [[codeblock]]
    // const Scalar minPressureGradient = getParam<Scalar>("Problem.MinimumPressureGradient", 1e1);
    // const Scalar maxPressureGradient = getParam<Scalar>("Problem.MaximumPressureGradient", 1e10);

    // if (!(minPressureGradient < maxPressureGradient))
    //     DUNE_THROW(Dune::InvalidStateException, "Maximum pressure gradient must be greater than minimum pressure gradient");

    // const int numberOfSamples = getParam<int>("Problem.NumberOfPressureGradients", 1);
    // // [[/codeblock]]

    // // Iterate over all directions specified before, apply several pressure gradients, calculate the mass flux
    // // and finally determine the the upscaled properties.
    // // [[codeblock]]
    // const auto directions = getParam<std::vector<std::size_t>>("Problem.Directions", std::vector<std::size_t>{0, 1, 2});
    // // upscalingHelper.setDirections(directions);

    // for (int dimIdx : directions)
    // {
    //     // set the direction in which the pressure gradient will be applied
    //     problem->setDirection(dimIdx);
    //     // ### The actual upscaling procedure
    //     // #### Instantiate the upscaling helper
    //     // [[codeblock]]

    //     upscalingHelper.setDirections(directions);
    //     // [[/codeblock]]

    //     // Set the side lengths used for applying the pressure gradient and calculating the REV outflow area.
    //     // One can either specify these values manually (usually more accurate) or let the UpscalingHelper struct
    //     // determine it automatically based on the network's bounding box.
    //     // [[codeblock]]
    //     const auto sideLengths = [&]()
    //     {
    //         using GlobalPosition = typename GridGeometry::GlobalCoordinate;
    //         if (hasParam("Problem.SideLength"))
    //             return getParam<GlobalPosition>("Problem.SideLength");
    //         else
    //             return upscalingHelper.getSideLengths(*gridGeometry);
    //     }();

    //     problem->setSideLengths(sideLengths);

    //     for (std::size_t step = 0; step<staticFlowProperties.size(); ++step)
    //     {
    //         problem->importTransmissibility(staticFlowProperties[step].throatTransmissibility);
    //         // for (int i = 0; i < 1; i++)
    //         // {
    //         int i = 0;
    //         // reset the solution
    //         x = 0;

    //         // set the pressure gradient to be applied
    //         Scalar pressureGradient = maxPressureGradient * std::exp(i + 1 - numberOfSamples);
    //         if (i == 0)
    //             pressureGradient = std::min(minPressureGradient, pressureGradient);

    //         problem->setPressureGradient(pressureGradient);

    //         // solve problem
    //         nonLinearSolver.solve(x);

    //         // set the sample points
    //         const Scalar totalFluidMassFlux = boundaryFlux.getFlux(std::vector<int>{ problem->outletPoreLabel() })[0];
    //         upscalingHelper.setDataPoints(*problem, totalFluidMassFlux, staticFlowProperties[step]);
    //         // }

    //         // write a vtu file for the given direction for the last sample
    //         vtkWriter.write(step);

    //     }

    //     // calculate and report the upscaled properties
    //     constexpr bool isCreepingFlow = 0;
    //     // upscalingHelper.calculateUpscaledProperties(*problem, isCreepingFlow);
    //     upscalingHelper.report(isCreepingFlow);

    //     // compare the Darcy permeability with reference data if provided in input file and report in case of inconsistency
    //     static const auto referenceData = getParam<std::vector<Scalar>>("Problem.ReferencePermeability", std::vector<Scalar>{});
    //     if (!referenceData.empty())
    //         upscalingHelper.compareWithReference(referenceData);

    //     upscalingHelper.writePlotDataToFile(dimIdx, phaseIdx);
    // }
};

template<class GridManager>
const auto runStaticProblem(GridManager& gridManager, int step)
{
    using TypeTag = Properties::TTag::PNMTWOPStatic;

    const auto& leafGridView = gridManager.grid().leafGridView();
    auto gridData = gridManager.getGridData();

    //////////////////////////////////////////////////////////////////////////////

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView, *gridData);


    using Problem = GetPropType<TypeTag, Properties::Problem>;
    Problem problem(*gridGeometry);

    Dune::VTKSequenceWriter<GridGeometry::GridView> sequenceWriter = std::move(problem.sequenceWriter());

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    PoreNetwork::TwoPStaticDrainage<GridGeometry, Scalar> drainageModel = std::move(problem.drainageModelStatic());


    problem.pcSwStatic(drainageModel, step);
    sequenceWriter.write(step);

    return problem.sw();
}

} // end namespace Dumux
// [[/codeblock]]

int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTagStatic = Properties::TTag::PNMTWOPStatic;
    // using TypeTagLiquid = Properties::TTag::PNMUpscalingCreepingFlowLiquid;
    // using TypeTagGas = Properties::TTag::PNMUpscalingCreepingFlowGas;

    using Scalar = GetPropType<TypeTagStatic, Properties::Scalar>;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    ////////////////////////////////////////////////////////////
    // parse the command line arguments and input file
    ////////////////////////////////////////////////////////////
    // parse command line arguments
    Parameters::init(argc, argv);

    using GridManager = PoreNetwork::GridManager<3>;
    GridManager gridManager;
    gridManager.init();


    auto numSteps = getParam<int>("Problem.NumSteps");
    // for (int step = 0; step < numSteps + 1; ++step)
    // {
        int step = 4;
        const auto sw = runStaticProblem(gridManager, step);
        runSinglePhaseUpscaling(gridManager, sw, step);
    // }

    // UpscalingHelper<Scalar>::plot();
    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////
    // timeLoop->finalize(leafGridView.comm());

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }
}
