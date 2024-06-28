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

template<class TypeTag, class GridManager, class StaticProperties>
void runSinglePhaseUpscaling(GridManager& gridManager, const StaticProperties& staticProperties)
{
    constexpr int phaseIdx = std::is_same_v<TypeTag, Properties::TTag::PNMUpscalingCreepingFlowLiquid> ? 0:1;
    const auto& staticFlowProperties = staticProperties;
    // ### Create the grid and the grid geometry
    // [[codeblock]]
    // The grid manager can be used to create a grid from the input file
    // using GridManager = PoreNetwork::GridManager<3>;
    // GridManager gridManager;
    // gridManager.init();

    // We compute on the leaf grid view.
    const auto& leafGridView = gridManager.grid().leafGridView();

    // instantiate the grid geometry
    auto gridData = gridManager.getGridData();
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView, *gridData);
    // [[/codeblock]]

    // ### Initialize the problem and grid variables
    // [[codeblock]]
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    problem->importTransmissibility(staticFlowProperties[0].throatTransmissibility);

    // instantiate and initialize the discrete and exact solution vectors
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs()); // zero-initializes

    // instantiate and initialize the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);
    // [[/codeblock]]

    // ### Instantiate the solver
    // We use the `NewtonSolver` class, which is instantiated on the basis
    // of an assembler and a linear solver. When the `solve` function of the
    // `NewtonSolver` is called, it uses the assembler and linear
    // solver classes to assemble and solve the non-linear system.
    // [[codeblock]]
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);


    using LinearSolver = Dumux::UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();

    using NewtonSolver = NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);
    // [[/codeblock]]

    // ### Initialize VTK output
    using VtkOutputFields = GetPropType<TypeTag, Properties::IOFields>;
    using VtkWriter = PoreNetwork::VtkOutputModule<GridVariables, GetPropType<TypeTag, Properties::FluxVariables>, SolutionVector>;
    auto name = phaseIdx ? "Upscaling_gas":"Upscaling_water";
    VtkWriter vtkWriter(*gridVariables, x, name);
    VtkOutputFields::initOutputModule(vtkWriter);
    // specify the field type explicitly since it may not be possible
    // to deduce this from the vector size in a pore network
    vtkWriter.addField(gridGeometry->poreVolume(), "poreVolume", Vtk::FieldType::vertex);
    vtkWriter.addField(gridGeometry->throatCrossSectionalArea(), "throatCrossSectionalArea", Vtk::FieldType::element);

    // ### Prepare the upscaling procedure.
    // Set up a helper class to determine the total mass flux leaving the network
    const auto boundaryFlux = PoreNetwork::BoundaryFlux(*gridVariables, assembler->localResidual(), x);

    // // ### The actual upscaling procedure
    // // #### Instantiate the upscaling helper
    // // [[codeblock]]
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using UpscalingHelper = UpscalingHelper<Scalar>;
    UpscalingHelper upscalingHelper;
    // // [[/codeblock]]

    // // Set the side lengths used for applying the pressure gradient and calculating the REV outflow area.
    // // One can either specify these values manually (usually more accurate) or let the UpscalingHelper struct
    // // determine it automatically based on the network's bounding box.
    // // [[codeblock]]
    // const auto sideLengths = [&]()
    // {
    //     using GlobalPosition = typename GridGeometry::GlobalCoordinate;
    //     if (hasParam("Problem.SideLength"))
    //         return getParam<GlobalPosition>("Problem.SideLength");
    //     else
    //         return upscalingHelper.getSideLengths(*gridGeometry);
    // }();

    // pass the side lengths to the problem
    // problem->setSideLengths(sideLengths);
    // [[/codeblock]]

    // Get the maximum and minimum pressure gradient and the population of sample points specified in the input file
    // [[codeblock]]
    const Scalar minPressureGradient = getParam<Scalar>("Problem.MinimumPressureGradient", 1e1);
    const Scalar maxPressureGradient = getParam<Scalar>("Problem.MaximumPressureGradient", 1e10);

    if (!(minPressureGradient < maxPressureGradient))
        DUNE_THROW(Dune::InvalidStateException, "Maximum pressure gradient must be greater than minimum pressure gradient");

    const int numberOfSamples = getParam<int>("Problem.NumberOfPressureGradients", 1);
    // [[/codeblock]]

    // Iterate over all directions specified before, apply several pressure gradients, calculate the mass flux
    // and finally determine the the upscaled properties.
    // [[codeblock]]
    const auto directions = getParam<std::vector<std::size_t>>("Problem.Directions", std::vector<std::size_t>{0, 1, 2});
    // upscalingHelper.setDirections(directions);

    for (int dimIdx : directions)
    {
        // set the direction in which the pressure gradient will be applied
        problem->setDirection(dimIdx);
        // ### The actual upscaling procedure
        // #### Instantiate the upscaling helper
        // [[codeblock]]

        upscalingHelper.setDirections(directions);
        // [[/codeblock]]

        // Set the side lengths used for applying the pressure gradient and calculating the REV outflow area.
        // One can either specify these values manually (usually more accurate) or let the UpscalingHelper struct
        // determine it automatically based on the network's bounding box.
        // [[codeblock]]
        const auto sideLengths = [&]()
        {
            using GlobalPosition = typename GridGeometry::GlobalCoordinate;
            if (hasParam("Problem.SideLength"))
                return getParam<GlobalPosition>("Problem.SideLength");
            else
                return upscalingHelper.getSideLengths(*gridGeometry);
        }();

        problem->setSideLengths(sideLengths);

        for (std::size_t step = 0; step<staticFlowProperties.size(); ++step)
        {
            problem->importTransmissibility(staticFlowProperties[step].throatTransmissibility);
            // for (int i = 0; i < 1; i++)
            // {
            int i = 0;
            // reset the solution
            x = 0;

            // set the pressure gradient to be applied
            Scalar pressureGradient = maxPressureGradient * std::exp(i + 1 - numberOfSamples);
            if (i == 0)
                pressureGradient = std::min(minPressureGradient, pressureGradient);

            problem->setPressureGradient(pressureGradient);

            // solve problem
            nonLinearSolver.solve(x);

            // set the sample points
            const Scalar totalFluidMassFlux = boundaryFlux.getFlux(std::vector<int>{ problem->outletPoreLabel() })[0];
            upscalingHelper.setDataPoints(*problem, totalFluidMassFlux, staticFlowProperties[step]);
            // }

            // write a vtu file for the given direction for the last sample
            vtkWriter.write(step);

        }

        // calculate and report the upscaled properties
        constexpr bool isCreepingFlow = 0;
        // upscalingHelper.calculateUpscaledProperties(*problem, isCreepingFlow);
        upscalingHelper.report(isCreepingFlow);

        // compare the Darcy permeability with reference data if provided in input file and report in case of inconsistency
        static const auto referenceData = getParam<std::vector<Scalar>>("Problem.ReferencePermeability", std::vector<Scalar>{});
        if (!referenceData.empty())
            upscalingHelper.compareWithReference(referenceData);

        upscalingHelper.writePlotDataToFile(dimIdx, phaseIdx);
    }
};

template<class GridManager>
const auto runStaticProblem(GridManager& gridManager)
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

    auto numSteps = problem.numberOfSteps();
    // const auto& flowProperties = problem.pcSwStatic(drainageModel);
    for (int step = 0; step < numSteps + 1; ++step)
    {
        problem.pcSwStatic(drainageModel, step);
        sequenceWriter.write(step);
    }

    return problem.staticProperties();
}

} // end namespace Dumux
// [[/codeblock]]

int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTagStatic = Properties::TTag::PNMTWOPStatic;
    using TypeTagLiquid = Properties::TTag::PNMUpscalingCreepingFlowLiquid;
    using TypeTagGas = Properties::TTag::PNMUpscalingCreepingFlowGas;

    using Scalar = GetPropType<TypeTagGas, Properties::Scalar>;

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

    const auto flowPropertiesStatic = runStaticProblem(gridManager);
    runSinglePhaseUpscaling<TypeTagLiquid>(gridManager, flowPropertiesStatic);
    runSinglePhaseUpscaling<TypeTagGas>(gridManager, flowPropertiesStatic);

    UpscalingHelper<Scalar>::plot();
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
