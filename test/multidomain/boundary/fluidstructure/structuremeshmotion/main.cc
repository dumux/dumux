// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#include <config.h>
#include <iostream>
#include <cmath>
#include <utility>

#include <dune/common/hybridutilities.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/elementsolution.hh>

#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/istlsolvers.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/fvassembler.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_alu.hh>
#include <dumux/io/grid/gridmanager_ug.hh>
#include <dumux/io/grid/gridmanager_sub.hh>

#include <dumux/experimental/timestepping/newmarkbeta.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);

    // initialize parameter tree
    Parameters::init(argc, argv);

    // type tags
    using StructureTypeTag = Properties::TTag::Structure;
    using MeshMotionTypeTag = Properties::TTag::MeshMotion;
    using CommonTypeTag = Properties::TTag::StructureMeshMotionCommon;

    using HostGrid = typename CommonTypeTag::HostGrid;
    GridManager<HostGrid> gridManager;
    gridManager.init();

    const auto& gridData = gridManager.getGridData();
    GridManager<GetPropType<CommonTypeTag, Properties::Grid>> gridManagerStructure, gridManagerFluid;
    gridManagerStructure.init(gridManager.grid(), [&](const auto& element) { return gridData->getElementDomainMarker(element) == 1; });
    gridManagerFluid.init(gridManager.grid(), [&](const auto& element){ return gridData->getElementDomainMarker(element) == 2; });

    // we compute on the leaf grid view
    const auto& leafGridViewStructure = gridManagerStructure.grid().leafGridView();
    const auto& leafGridViewFluid = gridManagerFluid.grid().leafGridView();

    // create the finite volume grid geometries and share the basic grid geometry
    using StructureGridGeometry = GetPropType<StructureTypeTag, Properties::GridGeometry>;
    using MeshMotionGridGeometry = GetPropType<MeshMotionTypeTag, Properties::GridGeometry>;
    auto structureGridGeometry = std::make_shared<StructureGridGeometry>(leafGridViewStructure);
    auto meshMotionGridGeometry = std::make_shared<MeshMotionGridGeometry>(leafGridViewFluid);

    // coupling manager coupling flow balance with the rest of the model
    // deals with the different discretization schemes
    using CouplingManager = GetPropType<StructureTypeTag, Properties::CouplingManager>;
    auto couplingManager = std::make_shared<CouplingManager>(structureGridGeometry, meshMotionGridGeometry);

    // the problem (initial and boundary conditions)
    using StructureProblem = GetPropType<StructureTypeTag, Properties::Problem>;
    auto structureProblem = std::make_shared<StructureProblem>(structureGridGeometry, couplingManager);

    using MeshMotionProblem = GetPropType<MeshMotionTypeTag, Properties::Problem>;
    auto meshMotionProblem = std::make_shared<MeshMotionProblem>(meshMotionGridGeometry, couplingManager);

    // the solution vector
    constexpr auto structureIdx = CouplingManager::structureIdx;
    constexpr auto meshMotionIdx = CouplingManager::meshMotionIdx;
    using Traits = MultiDomainTraits<StructureTypeTag, MeshMotionTypeTag>;
    using SolutionVector = typename Traits::SolutionVector;
    SolutionVector x;
    structureProblem->applyInitialSolution(x[structureIdx]);
    meshMotionProblem->applyInitialSolution(x[meshMotionIdx]);
    auto xOld = x;

    // the grid variables
    using StructureGridVariables = GetPropType<StructureTypeTag, Properties::GridVariables>;
    auto structureGridVariables = std::make_shared<StructureGridVariables>(structureProblem, structureGridGeometry);

    using MeshMotionGridVariables = GetPropType<MeshMotionTypeTag, Properties::GridVariables>;
    auto meshMotionGridVariables = std::make_shared<MeshMotionGridVariables>(meshMotionProblem, meshMotionGridGeometry);

    // initialize the coupling stencils
    couplingManager->init(structureProblem, meshMotionProblem, x);
    structureGridVariables->init(x[structureIdx]);
    meshMotionGridVariables->init(x[meshMotionIdx]);

    // VTK output
    using SolutionVectorStructure = std::decay_t<decltype(x[structureIdx])>;
    VtkOutputModule<StructureGridVariables, SolutionVectorStructure> vtkWriterStructure(*structureGridVariables, x[structureIdx], "structure");
    vtkWriterStructure.addVolumeVariable([](const auto& v){
        return Dune::FieldVector<double, 2>{v.displacement(0), v.displacement(1)};
    }, "d");
    std::vector<int> couplingMarkerStructure(structureGridGeometry->gridView().size(StructureGridGeometry::GridView::dimension), 0);
    for (std::size_t dofIdx = 0; dofIdx < couplingMarkerStructure.size(); ++dofIdx)
        couplingMarkerStructure[dofIdx] = couplingManager->isCoupledDof(dofIdx, CouplingManager::structureIdx) ? 1 : 0;
    vtkWriterStructure.addField(couplingMarkerStructure, "couplingMarker");
    vtkWriterStructure.write(0.0);

    using SolutionVectorMeshMotion = std::decay_t<decltype(x[meshMotionIdx])>;
    VtkOutputModule<MeshMotionGridVariables, SolutionVectorMeshMotion> vtkWriterFluid(*meshMotionGridVariables, x[meshMotionIdx], "fluid");
    vtkWriterFluid.addVolumeVariable([](const auto& v){
        return Dune::FieldVector<double, 2>{v.priVar(0), v.priVar(1)};
    }, "d");
    std::vector<int> couplingMarkerFluid(meshMotionGridGeometry->gridView().size(MeshMotionGridGeometry::GridView::dimension), 0);
    for (std::size_t dofIdx = 0; dofIdx < couplingMarkerFluid.size(); ++dofIdx)
        couplingMarkerFluid[dofIdx] = couplingManager->isCoupledDof(dofIdx, CouplingManager::meshMotionIdx) ? 1 : 0;
    vtkWriterFluid.addField(couplingMarkerFluid, "couplingMarker");
    vtkWriterFluid.write(0.0);

    // get some time loop parameters
    using Scalar = GetPropType<StructureTypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto dt = getParam<Scalar>("TimeLoop.DtInitial");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");

    // instantiate time loop
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    auto newmarkBeta = std::make_shared<Experimental::NewmarkBeta<Scalar, SolutionVectorStructure>>();
    newmarkBeta->initialize(x[structureIdx]);
    structureProblem->setNewmarkScheme(newmarkBeta);

    // the assembler with time loop for a transient problem
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(
        std::make_tuple(structureProblem, meshMotionProblem),
        std::make_tuple(structureGridGeometry, meshMotionGridGeometry),
        std::make_tuple(structureGridVariables, meshMotionGridVariables),
        couplingManager, timeLoop, xOld
    );

    // the linear solver
    using LAT = LinearAlgebraTraitsFromAssembler<Assembler>;
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LAT>;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = Dumux::MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    auto nonLinearSolver = std::make_shared<NewtonSolver>(assembler, linearSolver, couplingManager);

    const int vtkInterval = getParam<int>("VTKOutput.Every", 5);

    timeLoop->start(); do
    {
        // linearize & solve
        nonLinearSolver->solve(x);

        // update the solution in the time stepping scheme
        newmarkBeta->update(dt, x[structureIdx]);

        // make the new solution the old solution
        xOld = x;
        structureGridVariables->advanceTimeStep();
        meshMotionGridVariables->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // write VTK output
        if (timeLoop->timeStepIndex() % vtkInterval == 0 || timeLoop->finished())
        {
            vtkWriterStructure.write(timeLoop->time());
            vtkWriterFluid.write(timeLoop->time());
        }

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by the newton solver
        //timeLoop->setTimeStepSize(nonLinearSolver->suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridViewStructure.comm());

    return 0;
}
