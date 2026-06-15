// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BedloadTests
 * \brief Channel bend test for the bedload transport model.
 *
 * This test is based on run number four of the experiments presented in
 * Yen and Lee (1995) "Bed topography and sediment sorting in channel bend
 * with unsteady flow", Journal of Hydraulic Engineering.
 *
 * In addition to verifying the core functionality (as covered by the bump test),
 * this case tests further features of the bedload transport model, including:
 * - secondary current approach,
 * - the layer model,
 * - mixed-size bedload transport,
 * - the Meyer-Peter–Müller bedload transport formula.
 */
#include <config.h>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/io/grid/gridmanager_ug.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/istlsolverfactorybackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/assembly/fvassembler.hh>

#include "properties.hh"


int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tags for this problem
    using ShallowWaterTypeTag = Properties::TTag::ShallowWaterSub;
    using BedloadTypeTag = Properties::TTag::BedloadSub;

    // maybe initialize MPI and/or multithreading backend
    initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // create a grid (defined in params.input)
    // As the grids for the shallow water and the bedload model are identical it does not matter if the ShallowWaterTypeTag or the BedloadTypeTag is used
    GridManager<GetPropType<ShallowWaterTypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometries
    using ShallowWaterGridGeometry = GetPropType<ShallowWaterTypeTag, Properties::GridGeometry>;
    using BedloadGridGeometry = GetPropType<BedloadTypeTag, Properties::GridGeometry>;
    auto shallowWaterGridGeometry = std::make_shared<ShallowWaterGridGeometry>(leafGridView);
    auto bedloadGridGeometry = std::make_shared<BedloadGridGeometry>(leafGridView);

    // read the initial condition, which is the solution of a steady state simulation without bedload transport
    using Scalar = GetPropType<ShallowWaterTypeTag, Properties::Scalar>;
    const auto gridData = gridManager.getGridData();
    const auto& gridView = gridManager.grid().leafGridView();

    std::vector<Scalar> initialBedSurface(gridView.size(0));
    std::vector<Scalar> initialWaterDepth(gridView.size(0));
    std::vector<Scalar> initialVelocityX(gridView.size(0));
    std::vector<Scalar> initialVelocityY(gridView.size(0));
    std::vector<Scalar> curvatureRadius(gridView.size(0));
    for (const auto& element : elements(gridView))
    {
        const auto eIdx = gridView.indexSet().index(element);
        initialBedSurface[eIdx] = gridData->getParameter(element, "bedSurface");
        initialWaterDepth[eIdx] = gridData->getParameter(element, "h");
        initialVelocityX[eIdx] = gridData->getParameter(element, "u");
        initialVelocityY[eIdx] = gridData->getParameter(element, "v");
        curvatureRadius[eIdx] = gridData->getParameter(element, "curvatureRadius");
    }

    // the problems (initial and boundary conditions)
    using ShallowWaterProblem = GetPropType<ShallowWaterTypeTag, Properties::Problem>;
    auto shallowWaterSpatialParams = std::make_shared<typename ShallowWaterProblem::SpatialParams>(shallowWaterGridGeometry, initialBedSurface);
    auto shallowWaterProblem = std::make_shared<ShallowWaterProblem>(shallowWaterGridGeometry, shallowWaterSpatialParams, initialWaterDepth, initialVelocityX, initialVelocityY);
    using BedloadProblem = GetPropType<BedloadTypeTag, Properties::Problem>;
    auto bedloadSpatialParams = std::make_shared<typename BedloadProblem::SpatialParams>(bedloadGridGeometry, initialBedSurface, initialWaterDepth,
                                                                                         initialVelocityX, initialVelocityY, curvatureRadius);
    auto bedloadProblem = std::make_shared<BedloadProblem>(bedloadGridGeometry, bedloadSpatialParams, initialBedSurface);

    // the solution vectors
    using ShallowWaterSolutionVector = GetPropType<ShallowWaterTypeTag, Properties::SolutionVector>;
    using BedloadSolutionVector = GetPropType<BedloadTypeTag, Properties::SolutionVector>;
    ShallowWaterSolutionVector shallowWaterX(shallowWaterGridGeometry->numDofs());
    BedloadSolutionVector bedloadX(bedloadGridGeometry->numDofs());
    shallowWaterProblem->applyInitialSolution(shallowWaterX);
    bedloadProblem->applyInitialSolution(bedloadX);
    ShallowWaterSolutionVector shallowWaterXOld = shallowWaterX;
    BedloadSolutionVector bedloadXOld = bedloadX;

    // the grid variables
    using ShallowWaterGridVariables = GetPropType<ShallowWaterTypeTag, Properties::GridVariables>;
    using BedloadGridVariables = GetPropType<BedloadTypeTag, Properties::GridVariables>;
    auto shallowWaterGridVariables = std::make_shared<ShallowWaterGridVariables>(shallowWaterProblem, shallowWaterGridGeometry);
    auto bedloadGridVariables = std::make_shared<BedloadGridVariables>(bedloadProblem, bedloadGridGeometry);
    shallowWaterGridVariables->init(shallowWaterX);
    bedloadGridVariables->init(bedloadX);

    // get some time loop parameters
    using Scalar = GetPropType<ShallowWaterTypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto dtPlot = getParam<Scalar>("TimeLoop.DtPlot");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // initialize the vtk output module
    using IOFields = GetPropType<BedloadTypeTag, Properties::IOFields>;

    // write the initial state to vtk output
    VtkOutputModule<BedloadGridVariables, BedloadSolutionVector> vtkWriter(*bedloadGridVariables, bedloadX, bedloadProblem->name());
    IOFields::initOutputModule(vtkWriter);
    vtkWriter.write(0.0);

    // instantiate time loop
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // the assembler
    using ShallowWaterAssembler = FVAssembler<ShallowWaterTypeTag, DiffMethod::numeric>;
    using BedloadAssembler = FVAssembler<BedloadTypeTag, DiffMethod::numeric>;
    auto shallowWaterAssembler = std::make_shared<ShallowWaterAssembler>(shallowWaterProblem, shallowWaterGridGeometry, shallowWaterGridVariables, timeLoop, shallowWaterXOld);
    auto bedloadAssembler = std::make_shared<BedloadAssembler>(bedloadProblem, bedloadGridGeometry, bedloadGridVariables, timeLoop, bedloadXOld);

    // the linear solver (parallel)
    using ShallowWaterLinearSolver = IstlSolverFactoryBackend<LinearSolverTraits<ShallowWaterGridGeometry>, LinearAlgebraTraitsFromAssembler<ShallowWaterAssembler>>;
    using BedloadLinearSolver = IstlSolverFactoryBackend<LinearSolverTraits<BedloadGridGeometry>, LinearAlgebraTraitsFromAssembler<BedloadAssembler>>;
    auto shallowWaterLinearSolver = std::make_shared<ShallowWaterLinearSolver>(leafGridView, shallowWaterGridGeometry->dofMapper());
    auto bedloadLinearSolver = std::make_shared<BedloadLinearSolver>(leafGridView, bedloadGridGeometry->dofMapper());

    // the non-linear solver
    using ShallowWaterNewtonSolver = NewtonSolver<ShallowWaterAssembler, ShallowWaterLinearSolver>;
    using BedloadNewtonSolver = NewtonSolver<BedloadAssembler, BedloadLinearSolver>;
    ShallowWaterNewtonSolver shallowWaterNonLinearSolver(shallowWaterAssembler, shallowWaterLinearSolver);
    BedloadNewtonSolver bedloadNonLinearSolver(bedloadAssembler, bedloadLinearSolver);

    // set some check points
    Scalar const eps = 1e-8;
    for (Scalar i=dtPlot; i<tEnd + eps; i+=dtPlot)
    {
        timeLoop->setCheckPoint(i);
    }
    timeLoop->setCheckPoint(tEnd);

    // time loop
    timeLoop->start(); do
    {
        // update boundary conditions
        shallowWaterProblem->updateBoundaryConditions(timeLoop->time());

        // solve the shallow water system
        if (mpiHelper.rank() == 0) { std::puts("Solve the shallow water equation system:"); }
        shallowWaterNonLinearSolver.solve(shallowWaterX, *timeLoop);

        // update shallow water variables
        auto shallowWaterCoupledVariables = shallowWaterProblem->getCoupledVariables(shallowWaterX, *shallowWaterGridVariables);
        bedloadSpatialParams->updateCoupledVariables(shallowWaterCoupledVariables["h"], shallowWaterCoupledVariables["u"], shallowWaterCoupledVariables["v"],
                                                     shallowWaterCoupledVariables["bottomShearStressX"], shallowWaterCoupledVariables["bottomShearStressY"]);

        // solve the bedload system
        if (mpiHelper.rank() == 0) { std::puts("Solve the bedload equation system:"); }
        bedloadNonLinearSolver.solve(bedloadX,*timeLoop);

        // update bed suface
        auto bedloadCoupledVariables = bedloadProblem->getCoupledVariables(bedloadX, *bedloadGridVariables);
        shallowWaterSpatialParams->updateCoupledVariables(bedloadCoupledVariables["bedSurface"]);

        // update the layer model
        bedloadProblem->updateLayerModel(bedloadX, *bedloadGridVariables);

        // make the new solution the old solution
        shallowWaterXOld = shallowWaterX;
        bedloadXOld = bedloadX;
        shallowWaterGridVariables->advanceTimeStep();
        bedloadGridVariables->advanceTimeStep();

        // advance the time loop to the next step
        timeLoop->advanceTimeStep();

        // write vtk output
        if (timeLoop->isCheckPoint())
        {
            vtkWriter.write(timeLoop->time());
        }

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by newton controller
        using std::min;
        timeLoop->setTimeStepSize(min(shallowWaterNonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()),
                                      bedloadNonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize())));

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
    }

    return 0;
}
