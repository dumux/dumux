// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 *
 * \brief Two-phase relative-permeability upscaling with a pore network model.
 *
 * A quasi-static two-phase invasion-percolation sweep (TwoPStatic, driven by
 * PNMProblemStatic) determines, over a full capillary-pressure cycle (drainage from
 * Problem.InitialPc up to Problem.FinalPc, then imbibition back down), the invasion state and
 * phase conductance of every throat at each step. For each step, a single-phase Darcy
 * (CreepingFlow) solve is then run separately for the wetting and the non-wetting phase on that
 * frozen conductance field, and the resulting mass fluxes are converted into apparent phase
 * permeabilities, giving relative-permeability curves that capture drainage/imbibition
 * hysteresis.
 */
#include <config.h>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/discretization/porenetwork/gridgeometry.hh>
#include <dumux/io/grid/porenetwork/gridmanager.hh>
#include <dumux/porenetwork/common/boundaryflux.hh>
#include <dumux/porenetwork/common/pnmvtkoutputmodule.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/assembly/fvassembler.hh>

#include "upscalinghelper.hh"
#include "properties.hh"

namespace Dumux {

/*!
 * \brief Runs the quasi-static two-phase invasion-percolation sweep for the whole saturation
 *        range, caching the per-throat phase transmissibilities of every step.
 *
 * \param gridManager The pore-network grid manager
 * \param upscalingHelper Receives the network's total pore volume for the porosity report
 */
template<class GridManager>
auto runStaticProblem(GridManager& gridManager,
                      UpscalingHelper<GetPropType<Properties::TTag::PNMTWOPStatic, Properties::Scalar>>& upscalingHelper)
{
    using TypeTag = Properties::TTag::PNMTWOPStatic;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    const auto& leafGridView = gridManager.grid().leafGridView();
    auto gridData = gridManager.getGridData();
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView, *gridData);

    Problem problem(*gridGeometry);

    auto sequenceWriter = problem.sequenceWriter();

    for (int step = 0; step < problem.totalNumberOfSteps(); ++step)
    {
        problem.pcSwStatic(step);
        sequenceWriter.write(step);
    }

    if (getParam<bool>("Problem.PlotPcS", false))
        problem.plot();

    upscalingHelper.setTotalPoreVolume(problem.totalPoreVolume());
    upscalingHelper.setApexStep(problem.numberOfSteps());
    return problem.staticProperties();
}

/*!
 * \brief Runs the single-phase Darcy upscaling solve for one fluid (liquid or gas), looking
 *        up the throat conductances precomputed by the two-phase model at every saturation step.
 *
 * \tparam TypeTag Selects liquid or gas via its FluidSystem
 * \param gridManager The pore-network grid manager
 * \param staticFlowProperties The per-step results of the two-phase sweep (see runStaticProblem)
 * \param upscalingHelper Collects the resulting apparent-permeability-vs-saturation data
 */
template<class TypeTag, class GridManager, class StaticProperties>
void runSinglePhaseUpscaling(GridManager& gridManager,
                             const StaticProperties& staticFlowProperties,
                             UpscalingHelper<GetPropType<TypeTag, Properties::Scalar>>& upscalingHelper)
{
    constexpr int phaseIdx = std::is_same_v<TypeTag, Properties::TTag::PNMUpscalingCreepingFlowLiquid> ? 0 : 1;

    const auto& leafGridView = gridManager.grid().leafGridView();
    auto gridData = gridManager.getGridData();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView, *gridData);

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs()); // zero-initializes

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // The `NewtonSolver` is instantiated on the basis of an assembler and a linear solver; its
    // `solve()` function uses them to assemble and solve the (here: linear) system.
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    using LinearSolver = Dumux::UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();

    using NewtonSolver = NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    using VtkOutputFields = GetPropType<TypeTag, Properties::IOFields>;
    using VtkWriter = PoreNetwork::VtkOutputModule<GridVariables, GetPropType<TypeTag, Properties::FluxVariables>, SolutionVector>;
    const std::string name = phaseIdx == 0 ? "upscaling_liquid" : "upscaling_gas";
    VtkWriter vtkWriter(*gridVariables, x, name);
    VtkOutputFields::initOutputModule(vtkWriter);
    // specify the field type explicitly since it may not be possible
    // to deduce this from the vector size in a pore network
    vtkWriter.addField(gridGeometry->poreVolume(), "poreVolume", Vtk::FieldType::vertex);
    vtkWriter.addField(gridGeometry->throatCrossSectionalArea(), "throatCrossSectionalArea", Vtk::FieldType::element);

    // set up a helper to determine the total mass flux leaving the network
    const auto boundaryFlux = PoreNetwork::BoundaryFlux(*gridVariables, assembler->localResidual(), x);

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const Scalar minPressureGradient = getParam<Scalar>("Problem.MinimumPressureGradient", 1e1);
    const Scalar maxPressureGradient = getParam<Scalar>("Problem.MaximumPressureGradient", 1e10);
    if (!(minPressureGradient < maxPressureGradient))
        DUNE_THROW(Dune::InvalidStateException, "Maximum pressure gradient must be greater than minimum pressure gradient");

    // only the flow direction matters here (relative permeability is upscaled along the same
    // direction the two-phase sweep used to establish the saturation state)
    const auto directions = getParam<std::vector<std::size_t>>("Problem.Directions", std::vector<std::size_t>{0});
    for (int dimIdx : directions)
    {
        problem->setDirection(dimIdx);

        // Set the side lengths used for applying the pressure gradient and calculating the REV
        // outflow area. One can either specify these values manually (usually more accurate) or
        // let the UpscalingHelper struct determine them automatically from the network's
        // bounding box.
        const auto sideLengths = [&]()
        {
            using GlobalPosition = typename GridGeometry::GlobalCoordinate;
            if (hasParam("Problem.SideLength"))
                return getParam<GlobalPosition>("Problem.SideLength");
            else
                return upscalingHelper.getSideLengths(*gridGeometry);
        }();
        problem->setSideLengths(sideLengths);

        // start this phase's data for this direction fresh (the helper instance is shared
        // across both phases, see runSinglePhaseUpscaling's doc comment)
        upscalingHelper.resetDataPoints(dimIdx);

        for (std::size_t step = 0; step < staticFlowProperties.size(); ++step)
        {
            problem->importTransmissibility(staticFlowProperties[step].throatTransmissibility);
            x = 0.0;

            // apply the smallest available pressure gradient: for the creeping-flow regime used
            // here, the apparent permeability is independent of the pressure gradient, so a
            // single, small sample point per saturation step is sufficient
            problem->setPressureGradient(std::min(minPressureGradient, maxPressureGradient));

            nonLinearSolver.solve(x);

            const Scalar totalFluidMassFlux = boundaryFlux.getFlux(std::vector<int>{problem->outletPoreLabel()})[0];
            upscalingHelper.setDataPoints(*problem, totalFluidMassFlux, staticFlowProperties[step]);

            vtkWriter.write(step);
        }

        upscalingHelper.networkPorosity(*problem);
        upscalingHelper.writePlotDataToFile(dimIdx, phaseIdx, getParam<std::string>("Problem.Name"));
    }
}

} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTagLiquid = Properties::TTag::PNMUpscalingCreepingFlowLiquid;
    using TypeTagGas = Properties::TTag::PNMUpscalingCreepingFlowGas;
    using Scalar = GetPropType<TypeTagGas, Properties::Scalar>;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and the input file
    Parameters::init(argc, argv);

    using GridManager = PoreNetwork::GridManager<3>;
    GridManager gridManager;
    gridManager.init();

    UpscalingHelper<Scalar> upscalingHelper;
    const auto staticFlowProperties = runStaticProblem(gridManager, upscalingHelper);
    runSinglePhaseUpscaling<TypeTagLiquid>(gridManager, staticFlowProperties, upscalingHelper);
    runSinglePhaseUpscaling<TypeTagGas>(gridManager, staticFlowProperties, upscalingHelper);

    upscalingHelper.reportPorosity();
    if (getParam<bool>("Problem.PlotPcS", false))
        upscalingHelper.plot();

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
}
