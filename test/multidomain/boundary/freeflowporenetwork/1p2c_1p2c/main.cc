// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoundaryTests
 * \brief A test problem for the coupled Free-Flow/PNM problem (1p2c).
 */

#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>

#include <dumux/assembly/diffmethod.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/freeflow/navierstokes/velocityoutput.hh>
#include <dumux/io/grid/porenetwork/gridmanager.hh>
#include <dumux/io/vtk/intersectionwriter.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/porenetwork/common/pnmvtkoutputmodule.hh>
#include <dumux/multidomain/boundary/freeflowporenetwork/snappygridmanager.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include "properties.hh"

template<class GridGeometry, class GridVariables, class SolutionVector>
void updateVelocities(
    std::vector<Dune::FieldVector<double, 2>>& faceVelocity,
    const GridGeometry& gridGeometry,
    const GridVariables& gridVariables,
    const SolutionVector& x
){
    auto fvGeometry = localView(gridGeometry);
    auto elemVolVars = localView(gridVariables.curGridVolVars());
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        fvGeometry.bind(element);
        elemVolVars.bind(element, fvGeometry, x);

        for (const auto& scv : scvs(fvGeometry))
            faceVelocity[scv.dofIndex()][scv.dofAxis()] = elemVolVars[scv].velocity();
    }
}

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
    using FreeFlowMomentumTypeTag = Properties::TTag::FreeFlowOnePNCMomentum;
    using FreeFlowMassTypeTag = Properties::TTag::FreeFlowOnePNCMass;
    using PoreNetworkTypeTag = Properties::TTag::PNMOnePNCModel;

    using PNMGridManager = Dumux::PoreNetwork::GridManager<2>;
    PNMGridManager pnmGridManager;
    pnmGridManager.init("PNM");

    using FreeFlowGridManager = Dumux::PoreNetwork::SnappyGridManager<2, PNMGridManager>;
    FreeFlowGridManager freeflowGridManager;
    freeflowGridManager.init(pnmGridManager.grid(), *(pnmGridManager.getGridData()), "FreeFlow");

    // create the free flow grid geometries
    using FreeFlowMomentumGridGeometry = GetPropType<FreeFlowMomentumTypeTag, Properties::GridGeometry>;
    auto freeFlowMomentumGridGeometry = std::make_shared<FreeFlowMomentumGridGeometry>(freeflowGridManager.grid().leafGridView());

    using FreeFlowMassGridGeometry = GetPropType<FreeFlowMassTypeTag, Properties::GridGeometry>;
    auto freeFlowMassGridGeometry = std::make_shared<FreeFlowMassGridGeometry>(freeflowGridManager.grid().leafGridView());

    // create the PNM grid geometry
    using PoreNetworkGridGeometry = GetPropType<PoreNetworkTypeTag, Properties::GridGeometry>;
    auto pnmGridData = pnmGridManager.getGridData();
    auto pnmGridGeometry = std::make_shared<PoreNetworkGridGeometry>(pnmGridManager.grid().leafGridView(), *pnmGridData);

    using CouplingManager = GetPropType<FreeFlowMomentumTypeTag, Properties::CouplingManager>;
    auto couplingManager = std::make_shared<CouplingManager>();

    // the problem (initial and boundary conditions)
    using FreeFlowMomentumProblem = GetPropType<FreeFlowMomentumTypeTag, Properties::Problem>;
    auto freeFlowMomentumProblem = std::make_shared<FreeFlowMomentumProblem>(freeFlowMomentumGridGeometry, couplingManager);
    using FreeFlowMassProblem = GetPropType<FreeFlowMassTypeTag, Properties::Problem>;
    auto freeFlowMassProblem = std::make_shared<FreeFlowMassProblem>(freeFlowMassGridGeometry, couplingManager);

    using PoreNetworkSpatialParams = GetPropType<PoreNetworkTypeTag, Properties::SpatialParams>;
    auto pnmSpatialParams = std::make_shared<PoreNetworkSpatialParams>(pnmGridGeometry);
    using PNMProblem = GetPropType<PoreNetworkTypeTag, Properties::Problem>;
    auto pnmProblem = std::make_shared<PNMProblem>(pnmGridGeometry, pnmSpatialParams, couplingManager);

    // the indices
    constexpr auto freeFlowMomentumIndex = CouplingManager::freeFlowMomentumIndex;
    constexpr auto freeFlowMassIndex = CouplingManager::freeFlowMassIndex;
    constexpr auto poreNetworkIndex = CouplingManager::poreNetworkIndex;

    // the solution vector
    using Traits = MultiDomainTraits<FreeFlowMomentumTypeTag, FreeFlowMassTypeTag, PoreNetworkTypeTag>;
    Traits::SolutionVector sol;
    sol[freeFlowMomentumIndex].resize(freeFlowMomentumGridGeometry->numDofs());
    sol[freeFlowMassIndex].resize(freeFlowMassGridGeometry->numDofs());
    sol[poreNetworkIndex].resize(pnmGridGeometry->numDofs());
    freeFlowMassProblem->applyInitialSolution(sol[freeFlowMassIndex]);
    freeFlowMomentumProblem->applyInitialSolution(sol[freeFlowMomentumIndex]);
    pnmProblem->applyInitialSolution(sol[poreNetworkIndex]);
    auto solOld = sol;

    // the grid variables
    using FreeFlowMomentumGridVariables = GetPropType<FreeFlowMomentumTypeTag, Properties::GridVariables>;
    auto freeFlowMomentumGridVariables = std::make_shared<FreeFlowMomentumGridVariables>(freeFlowMomentumProblem, freeFlowMomentumGridGeometry);
    using FreeFlowMassGridVariables = GetPropType<FreeFlowMassTypeTag, Properties::GridVariables>;
    auto freeFlowMassGridVariables = std::make_shared<FreeFlowMassGridVariables>(freeFlowMassProblem, freeFlowMassGridGeometry);

    using PoreNetworkgridVariables = GetPropType<PoreNetworkTypeTag, Properties::GridVariables>;
    auto pnmGridVariables = std::make_shared<PoreNetworkgridVariables>(pnmProblem, pnmGridGeometry);

    couplingManager->init(freeFlowMomentumProblem, freeFlowMassProblem, pnmProblem,
                          std::make_tuple(freeFlowMomentumGridVariables, freeFlowMassGridVariables, pnmGridVariables),
                          sol, solOld);

    freeFlowMomentumGridVariables->init(sol[freeFlowMomentumIndex]);
    freeFlowMassGridVariables->init(sol[freeFlowMassIndex]);
    pnmGridVariables->init(sol[poreNetworkIndex]);

    // We get some time loop parameters from the input file
    // and instantiate the time loop
    using Scalar = typename Traits::Scalar;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    const auto dt = getParam<Scalar>("TimeLoop.DtInitial");
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // the assembler for a stationary problem
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(freeFlowMomentumProblem, freeFlowMassProblem, pnmProblem),
                                                 std::make_tuple(freeFlowMomentumGridGeometry,
                                                                 freeFlowMassGridGeometry,
                                                                 pnmGridGeometry),
                                                 std::make_tuple(freeFlowMomentumGridVariables,
                                                                 freeFlowMassGridVariables,
                                                                 pnmGridVariables),
                                                 couplingManager, timeLoop, solOld);

    // the linear solver
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    // initialize the vtk output module
    VtkOutputModule vtkWriterFFMass(*freeFlowMassGridVariables, sol[freeFlowMassIndex], freeFlowMassProblem->name());
    GetPropType<FreeFlowMassTypeTag, Properties::IOFields>::initOutputModule(vtkWriterFFMass);
    vtkWriterFFMass.addVelocityOutput(std::make_shared<NavierStokesVelocityOutput<FreeFlowMassGridVariables>>());

    using PoreNetworkOutputModule = PoreNetwork::VtkOutputModule<PoreNetworkgridVariables, GetPropType<PoreNetworkTypeTag, Properties::FluxVariables>, std::decay_t<decltype(sol[poreNetworkIndex])>>;
    PoreNetworkOutputModule pmVtkWriter(*pnmGridVariables, sol[poreNetworkIndex], pnmProblem->name());
    GetPropType<PoreNetworkTypeTag, Properties::IOFields>::initOutputModule(pmVtkWriter);

    // write vtk output
    vtkWriterFFMass.write(0.0);
    pmVtkWriter.write(0.0);

    // timeloop
    timeLoop->start(); do
    {
        // solve the non-linear system with time step control
        nonLinearSolver.solve(sol, *timeLoop);

        // make the new solution the old solution.
        solOld = sol;
        freeFlowMassGridVariables->advanceTimeStep();
        freeFlowMomentumGridVariables->advanceTimeStep();
        pnmGridVariables->advanceTimeStep();

        // advance to the time loop to the next step.
        timeLoop->advanceTimeStep();

        // write vtk output for each time step.
        vtkWriterFFMass.write(timeLoop->time());
        pmVtkWriter.write(timeLoop->time());

        // report statistics of this time step.
        timeLoop->reportTimeStep();

        // set a new dt as suggested by the newton solver for the next time step.
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());

    std::vector<Dune::FieldVector<double, 2>> faceVelocity(sol[freeFlowMomentumIndex].size());
    updateVelocities(faceVelocity, *freeFlowMomentumGridGeometry, *freeFlowMomentumGridVariables, sol[freeFlowMomentumIndex]);

    ConformingIntersectionWriter faceVtk(freeFlowMomentumGridGeometry->gridView());
    faceVtk.addField(faceVelocity, "velocityVector");
    faceVtk.write("face_velocities");

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
