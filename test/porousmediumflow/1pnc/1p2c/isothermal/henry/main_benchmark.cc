// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePNCTests
 * \brief Adaptive, single-rank benchmark for the Henry problem of Fahs et al. (2016, WRR,
 *        doi:10.1002/2016WR019288): same problem/properties as main.cc, but on ALUGrid
 *        (via properties.hh's HenryFahsBenchmarkTest/HenryFahsCase2BenchmarkTest type
 *        tags) with h-adaptive refinement/coarsening (see adaptive/gridadaptindicator.hh)
 *        instead of a fixed grid.
 *
 * The saltwater/freshwater mixing front is tracked by HenryGridAdaptIndicator (max jump
 * in X^solute between neighboring elements, normalized by the field's global range --
 * see [Adaptive] in params_benchmark(_case2).input for the refine/coarsen tolerances and
 * MinLevel/MaxLevel). There is no phase mass/state to conserve here (OnePNC's X^solute is
 * a plain primary variable, no primary-variable switching), so griddatatransfer.hh is
 * correspondingly simpler than the two-phase reference it's modeled on
 * (dumux/porousmediumflow/2p/griddatatransfer.hh).
 *
 *   ./test_1p2c_henry_case1_benchmark params_benchmark_case1.input -Problem.Name run1
 */

#include <config.h>

#include <cstddef>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/timeloop.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/assembly/fvassembler.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_alu.hh>

#include <dumux/adaptive/adapt.hh>
#include <dumux/adaptive/markelements.hh>

#include "properties.hh"
#include "adaptive/gridadaptindicator.hh"
#include "adaptive/griddatatransfer.hh"

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

    // initialize parameter tree
    Parameters::init(argc, argv);

    // try to create a grid (from the given grid file or the input file)
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    // the problem (initial and boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x;
    problem->applyInitialSolution(x);
    auto xOld = x;

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // --- adaptive refinement indicator + h-adaptive data transfer ---
    // primary-variable index of the solute mass fraction (problem.hh already indexes
    // PrimaryVariables/NumEqVector with this directly)
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    static constexpr int soluteIdx = FluidSystem::soluteIdx;

    HenryGridAdaptIndicator<TypeTag> indicator(gridGeometry, soluteIdx);
    HenryFahsBoxGridDataTransfer<TypeTag> dataTransfer(gridGeometry, x, xOld);

    const auto refineTol = getParam<double>("Adaptive.RefineTolerance", 0.05);
    const auto coarsenTol = getParam<double>("Adaptive.CoarsenTolerance", 0.001);

    // Initial concentration is uniform (domain filled with seawater, see problem.hh's
    // initialAtPos()), so the indicator has nothing to refine on before the front develops.
    // Still calculate once so indicator.values() is correctly sized for the VTK field
    // registration below.
    indicator.calculate(x, refineTol, coarsenTol);

    // initialize the vtk output module
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields

    // "adaptIndicatorValue": indicator.values(), the max jump in X^solute across an
    // element's faces that operator() compares against refineBound()/coarsenBound().
    // "adaptMark": operator()'s decision for each element (1 refine / -1 coarsen / 0 keep).
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    std::vector<Scalar> adaptMarkField;
    const auto updateAdaptMarkField = [&]()
    {
        const auto& gv = gridGeometry->gridView();
        adaptMarkField.assign(gv.size(0), 0.0);
        for (const auto& element : elements(gv))
            adaptMarkField[gridGeometry->elementMapper().index(element)] = static_cast<Scalar>(indicator(element));
    };
    updateAdaptMarkField();
    vtkWriter.addField(indicator.values(), "adaptIndicatorValue", Vtk::FieldType::element);
    vtkWriter.addField(adaptMarkField, "adaptMark", Vtk::FieldType::element);

    vtkWriter.write(0.0);

    // Variable time steps: a small DtInitial (see params_benchmark_case1.input) keeps the
    // first, steepest-front Newton solves easier. Newton's suggestTimeStepSize() grows
    // Delta t back up (capped at MaxTimeStepSize) as convergence gets easier. Output is
    // still pinned to the MaxTimeStepSize cadence (see the periodic check point below),
    // so VTU output times land on a fixed time grid regardless of the internal, adaptive
    // Delta t -- no interpolation needed downstream (e.g. in post_processing.py) to
    // compare animations.
    const auto tEnd = getParam<double>("TimeLoop.TEnd");
    const auto dtInit = getParam<double>("TimeLoop.DtInitial");
    const auto maxDt = getParam<double>("TimeLoop.MaxTimeStepSize");
    auto timeLoop = std::make_shared<CheckPointTimeLoop<double>>(0.0, dtInit, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);
    timeLoop->setPeriodicCheckPoint(maxDt);

    // the assembler with time loop for instationary problem
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);

    // UMFPack: a direct, sequential solver. The point of this benchmark is checking
    // whether h-adaptive refinement/coarsening alone converges to the right solution.
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    NewtonSolver<Assembler, LinearSolver> nonLinearSolver(assembler, linearSolver);

    // wall-clock timing
    Dune::Timer timer;

    // time loop
    std::size_t stepIdx = 0;
    timeLoop->start(); do
    {
        // Adapt only between timesteps: at this point x == xOld (both hold the solution
        // at the end of the previous accepted timestep, see the "xOld = x" below), so
        // interpolating both through dataTransfer is consistent. Skipped on the very
        // first iteration (nothing to adapt to yet: the initial concentration is
        // uniform, see the indicator.calculate() call above).
        if (stepIdx > 0)
        {
            indicator.calculate(x, refineTol, coarsenTol);

            if (markElements(gridManager.grid(), indicator) && adapt(gridManager.grid(), dataTransfer))
            {
                gridVariables->updateAfterGridAdaption(x);
                assembler->updateAfterGridAdaption();

                // indicator.values() was sized/indexed for the pre-adapt mesh; element
                // indices shift on adapt, not just their count -- refresh before it's
                // used for adaptMarkField/VTK output below.
                indicator.calculate(x, refineTol, coarsenTol);
            }
        }

        // solve the non-linear system for this time step
        nonLinearSolver.solve(x, *timeLoop);

        // make the new solution the old solution
        xOld = x;
        gridVariables->advanceTimeStep();

        // advance the time loop to the next, adaptively sized step
        timeLoop->advanceTimeStep();

        // write vtk output only at the periodic check points set up above (every
        // MaxTimeStepSize seconds of simulated time), not every internal, variable step
        if (timeLoop->isCheckPoint() || timeLoop->finished())
        {
            updateAdaptMarkField();
            vtkWriter.write(timeLoop->time());
        }

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // let Newton suggest the next time step size (grows/shrinks with convergence
        // ease, capped at MaxTimeStepSize; setPeriodicCheckPoint above still shrinks it
        // further as needed to hit the next output time exactly)
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));
        ++stepIdx;

    } while (!timeLoop->finished());

    timeLoop->finalize(gridGeometry->gridView().comm());

    const auto elapsed = timer.elapsed();
    if (mpiHelper.rank() == 0)
        std::cout << "\n[benchmark] solver = UMFPack, total wall-clock time = " << elapsed << " s\n" << std::endl;

    // print dumux end message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/false);

    return 0;
}
