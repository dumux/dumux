// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
// ## The main program (`main.cc`)
// This file contains the main program flow. In this example, we use a single-phase
// Pore-Network-Model to evaluate the upscaled Darcy permeability of a given network.
// [[content]]
// ### Includes
// [[details]] includes
// [[codeblock]]
#include <config.h>

#include <iostream>

#include <algorithm>

#include <dune/common/float_cmp.hh> // for floating point comparison

#include <dumux/common/properties.hh> // for GetPropType
#include <dumux/common/parameters.hh> // for getParam
#include <dumux/common/initialize.hh>

#include <dumux/linear/seqsolverbackend.hh> // for ILU0BiCGSTABBackend
#include <dumux/linear/pdesolver.hh>        // for LinearPDESolver
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/assembly/fvassembler.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/porenetwork/common/pnmvtkoutputmodule.hh>

#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/grid/porenetwork/gridmanager.hh> // for pore-network grid
#include <dumux/porenetwork/common/boundaryflux.hh> // for getting the total mass flux leaving the network

#include "upscalinghelper.hh"
#include "properties.hh"
// [[/codeblock]]
// [[/details]]
//
// ### The driver function
//
// It depends on the template argument `TypeTag` if we run the example assuming
// a creeping flow regime or not. This is decided with the parameter
// `Problem.AssumeCreepingFlow` in the input file.
namespace Dumux {

template<class TypeTag>
void runExample()
{
    // ### Create the grid and the grid geometry
    // [[codeblock]]
    // The grid manager can be used to create a grid from the input file
    using GridManager = PoreNetwork::GridManager<3>;
    GridManager gridManager;
    gridManager.init();

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

    using LinearSolver = UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    using NewtonSolver = NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);
    // [[/codeblock]]

    // ### Initialize VTK output
    using VtkOutputFields = GetPropType<TypeTag, Properties::IOFields>;
    using VtkWriter = PoreNetwork::VtkOutputModule<GridVariables, GetPropType<TypeTag, Properties::FluxVariables>, SolutionVector>;
    VtkWriter vtkWriter(*gridVariables, x, problem->name());
    VtkOutputFields::initOutputModule(vtkWriter);
    // specify the field type explicitly since it may not be possible
    // to deduce this from the vector size in a pore network
    vtkWriter.addField(gridGeometry->poreVolume(), "poreVolume", VtkWriter::FieldType::vertex);
    vtkWriter.addField(gridGeometry->throatShapeFactor(), "throatShapeFactor", VtkWriter::FieldType::element);
    vtkWriter.addField(gridGeometry->throatCrossSectionalArea(), "throatCrossSectionalArea", VtkWriter::FieldType::element);

    // ### Prepare the upscaling procedure.
    // Set up a helper class to determine the total mass flux leaving the network
    const auto boundaryFlux = PoreNetwork::BoundaryFlux(*gridVariables, assembler->localResidual(), x);

    // ### The actual upscaling procedure
    // #### Instantiate the upscaling helper
    // [[codeblock]]
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using UpscalingHelper = UpscalingHelper<Scalar>;
    UpscalingHelper upscalingHelper;
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

    // pass the side lengths to the problem
    problem->setSideLengths(sideLengths);
    // [[/codeblock]]

    // Get the maximum and minimum pressure gradient and the population of sample points specified in the input file
    // [[codeblock]]
    const Scalar minPressureGradient = getParam<Scalar>("Problem.MinimumPressureGradient", 1e1);
    const Scalar maxPressureGradient = getParam<Scalar>("Problem.MaximumPressureGradient", 1e10);

    if (!(minPressureGradient < maxPressureGradient))
        throw std::runtime_error("maximum pressure gradient must be greater than minimum pressure gradient");

    const int numberOfSamples = getParam<int>("Problem.NumberOfPressureGradients", 1);
    // [[/codeblock]]

    // Iterate over all directions specified before, apply several pressure gradient, calculated the mass flux
    // and finally determine the the upscaled properties.
    // [[codeblock]]
    const auto directions = getParam<std::vector<std::size_t>>("Problem.Directions", std::vector<std::size_t>{0, 1, 2});
    upscalingHelper.setDirections(directions);
    for (int dimIdx : directions)
    {
        // set the direction in which the pressure gradient will be applied
        problem->setDirection(dimIdx);

        for (int i = 0; i < numberOfSamples; i++)
        {
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
            upscalingHelper.setDataPoints(*problem, totalFluidMassFlux);
        }

        // write a vtu file for the given direction for the last sample
        vtkWriter.write(dimIdx);
    }

    // calculate and report the upscaled properties
    constexpr bool isCreepingFlow = std::is_same_v<TypeTag, Properties::TTag::PNMUpscalingCreepingFlow>;
    upscalingHelper.calculateUpscaledProperties(*problem, isCreepingFlow);
    upscalingHelper.report(isCreepingFlow);

    // compare the Darcy permeability with reference data if provided in input file and report in case of inconsistency
    static const auto referenceData = getParam<std::vector<Scalar>>("Problem.ReferencePermeability", std::vector<Scalar>{});
    if (!referenceData.empty())
        upscalingHelper.compareWithReference(referenceData);

    // plot the results just for non-creeping flow
    // creeping flow would just result in a straight line (permeability is independent of the pressure gradient)
    if (!isCreepingFlow)
        upscalingHelper.plot();
};

} // end namespace Dumux
// [[/codeblock]]

// ### The main function
// [[details]] main
// [[codeblock]]
int main(int argc, char** argv)
{
    using namespace Dumux;

    // We parse the command line arguments.
    Parameters::init(argc, argv);

    // Convenience alias for the type tag of the problem.
    using CreepingFlowTypeTag = Properties::TTag::PNMUpscalingCreepingFlow;
    using NonCreepingFlowTypeTag = Properties::TTag::PNMUpscalingNonCreepingFlow;
    // // [[/codeblock]]

    // user decides whether creeping flow or non-creeping flow should be run
    if (getParam<bool>("Problem.AssumeCreepingFlow", false))
        runExample<CreepingFlowTypeTag>();
    else
        runExample<NonCreepingFlowTypeTag>();

    // program end, return with 0 exit code (success)
    return 0;
}
// [[/codeblock]]
// [[/details]]
// [[/content]]
