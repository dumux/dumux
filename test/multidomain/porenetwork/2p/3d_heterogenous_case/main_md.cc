// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 *
 * \brief Test for the two-phase pore-network model
 */
#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/io/grid/porenetwork/gridmanager.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/porenetwork/common/pnmvtkoutputmodule.hh>

#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/porenetwork/common/outletpcgradient.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    using PNMTypeTag = Properties::TTag::DrainageProblem;
    using ConstraintTypeTag = Properties::TTag::ConstraintProblem;

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

    //////////////////////////////////////////////////////////////////////
    // try to create a grid (from the given grid file or the input file)
    /////////////////////////////////////////////////////////////////////

    using GridManager = Dumux::PoreNetwork::GridManager<3>;
    GridManager gridManager;
    gridManager.init();

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();
    auto gridData = gridManager.getGridData();

    // create the finite volume grid geometry
    using PNMGridGeometry = GetPropType<PNMTypeTag, Properties::GridGeometry>;
    auto pnmGridGeometry = std::make_shared<PNMGridGeometry>(leafGridView, *gridData);
    using ConstraintGridGeometry = GetPropType<ConstraintTypeTag, Properties::GridGeometry>;
    auto constraintGridGeometry = std::make_shared<ConstraintGridGeometry>(leafGridView);

    // the coupling manager
    using Traits = MultiDomainTraits<PNMTypeTag, ConstraintTypeTag>;
    using CouplingManager = PNMConstraintCouplingManager< Traits >;
    auto couplingManager = std::make_shared<CouplingManager>();

    // the spatial parameters
    using PNMSpatialParams = GetPropType<PNMTypeTag, Properties::SpatialParams>;
    auto pnmSpatialParams = std::make_shared<PNMSpatialParams>(pnmGridGeometry);

    // the problem (boundary conditions)
    using PNMProblem = GetPropType<PNMTypeTag, Properties::Problem>;
    auto pnmProblem = std::make_shared<PNMProblem>(pnmGridGeometry, pnmSpatialParams, "Pnm", couplingManager);
    using ConstraintProblem = GetPropType<ConstraintTypeTag, Properties::Problem>;
    auto constraintProblem = std::make_shared<ConstraintProblem>(constraintGridGeometry, "Constraint", couplingManager);

    // the solution vectors
    using SolutionVector = typename Traits::SolutionVector;
    SolutionVector x;
    static const auto pnmId = Traits::template SubDomain<0>::Index();
    static const auto constraintId = Traits::template SubDomain<1>::Index();
    x[pnmId].resize(pnmGridGeometry->numDofs());
    x[constraintId].resize(constraintGridGeometry->numDofs());

    // initialize the coupling manager
    couplingManager->init(pnmProblem, constraintProblem, x);
    pnmProblem->applyInitialSolution(x[pnmId]);
    pnmProblem->calculateSumInletVolume();
    constraintProblem->applyInitialSolution(x[constraintId]);
    auto xOld = x;

    // the grid variables
    using PNMGridVariables = GetPropType<PNMTypeTag, Properties::GridVariables>;
    using ConstraintGridVariables = GetPropType<ConstraintTypeTag, Properties::GridVariables>;
    auto pnmGridVariables = std::make_shared<PNMGridVariables>(pnmProblem, pnmGridGeometry);
    auto constraintGridVariables = std::make_shared<ConstraintGridVariables>(constraintProblem, constraintGridGeometry);
    pnmGridVariables->init(x[pnmId]);
    constraintGridVariables->init(x[constraintId]);

    couplingManager->setGridVariables(std::make_tuple(pnmGridVariables, constraintGridVariables));

    // get some time loop parameters
    const auto dt = getParam<double>("TimeLoop.DtInitial");
    const auto checkPoints = getParam<std::vector<double>>("TimeLoop.CheckPoints");
    const auto tEnd = getParam<double>("TimeLoop.TEnd");

    // initialize the vtk output module
    using IOFieldsPNM = GetPropType<PNMTypeTag, Properties::IOFields>;
    PoreNetwork::VtkOutputModule<PNMGridVariables, GetPropType<PNMTypeTag, Properties::FluxVariables>, GetPropType<PNMTypeTag, Properties::SolutionVector>> pnmVtkWriter(*pnmGridVariables, x[pnmId], pnmProblem->name());
    IOFieldsPNM::initOutputModule(pnmVtkWriter); //! Add model specific output fields

    using IOFieldsConstraint = GetPropType<ConstraintTypeTag, Properties::IOFields>;
    Dumux::VtkOutputModule<ConstraintGridVariables, GetPropType<ConstraintTypeTag, Properties::SolutionVector>> constraintVtkWriter(*constraintGridVariables, x[constraintId], constraintProblem->name());
    IOFieldsConstraint::initOutputModule(constraintVtkWriter); //! Add model specific output fields

    pnmVtkWriter.write(0.0);
    constraintVtkWriter.write(0.0);

    // use zero pc gradient BC for this test case
    const auto outletCapPressureGradient = std::make_shared<Dumux::PoreNetwork::OutletCapPressureGradient<PNMGridVariables, GetPropType<PNMTypeTag, Properties::SolutionVector>>>(*pnmGridVariables, x[pnmId]);
    pnmProblem->outletCapPressureGradient(outletCapPressureGradient);

    // instantiate time loop
    auto timeLoop = std::make_shared<CheckPointTimeLoop<double>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(getParam<double>("TimeLoop.MaxTimeStepSize"));
    timeLoop->setCheckPoint(checkPoints.begin(), checkPoints.end());

    // the assembler
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric, /*implicit?*/true>;
    auto assembler = std::make_shared<Assembler>( std::make_tuple(pnmProblem, constraintProblem),
                                                  std::make_tuple(pnmGridGeometry, constraintGridGeometry),
                                                  std::make_tuple(pnmGridVariables, constraintGridVariables),
                                                  couplingManager,
                                                  timeLoop,
                                                  xOld);

    // the linear solver
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver newtonSolver(assembler, linearSolver, couplingManager);

    // The following update is needed because the constraint is formualted for each time step and needs to be updated afterwards
    // One could also change the constraint by accounting for previous time step data
    auto updateState = [&]()
    {
        const auto& state = pnmGridVariables->gridFluxVarsCache().invasionState();
        for (const auto& element : elements(pnmGridGeometry->gridView()))
        {
            const auto eIdx = constraintGridGeometry->elementMapper().index(element);
            x[constraintId][eIdx] = (double) state.invaded(element);
        }
    };

    // time loop
    timeLoop->start(); do
    {
        if(pnmGridVariables->gridFluxVarsCache().invasionState().hasChanged())
        {
            updateState();
            constraintGridVariables->update(x[constraintId]);
        }

        // try solving the non-linear system
        newtonSolver.solve(x, *timeLoop);

        // make the new solution the old solution
        xOld = x;
        pnmGridVariables->advanceTimeStep();
        constraintGridVariables->advanceTimeStep();

        if (pnmGridVariables->gridFluxVarsCache().invasionState().updateAfterTimeStep())
            pnmGridVariables->gridFluxVarsCache().invasionState().update(x[pnmId], pnmGridVariables->curGridVolVars(), pnmGridVariables->gridFluxVarsCache());


        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // write vtk output
        if(pnmProblem->shouldWriteOutput(timeLoop->timeStepIndex(), *pnmGridVariables))
        {
            pnmVtkWriter.write(timeLoop->time());
            constraintVtkWriter.write(timeLoop->time());
        }

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by newton solver
        timeLoop->setTimeStepSize(newtonSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());

    newtonSolver.report();
    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////
    timeLoop->finalize(leafGridView.comm());

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }
}
