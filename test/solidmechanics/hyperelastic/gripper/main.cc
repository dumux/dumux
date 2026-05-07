// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Gripper hydrogel actuator — time-dependent multidomain driver.
 *
 * Two-subdomain problem:
 *  - Subdomain 0 (momentum, PQ1BubbleHybrid): displacement \f$\mathbf{u}\f$.
 *  - Subdomain 1 (mass, Box/P1): pore pressure \f$p\f$.
 *
 * The monolithic system is assembled with `MultiDomainFVAssembler` and solved
 * with `MultiDomainNewtonSolver` (using UMFPack as linear backend).
 *
 * After each converged time step:
 *  1. `couplingManager->savePrevTimeStepJacobian()` stores J^n for the
 *     mass storage term.
 *  2. The time is forwarded to the mass problem for the swelling BC.
 */
#include <config.h>
#include <iostream>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/timeloop.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/istlsolvers.hh>

#include <dumux/experimental/timestepping/newmarkbeta.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_foam.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv);

    using MomentumTypeTag = Properties::TTag::GripperMomentum;
    using MassTypeTag = Properties::TTag::GripperMass;
    using MDTraits = GripperMDTraits; // defined in properties.hh

    using Grid = GetPropType<MomentumTypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();

    using MomentumGridGeometry = GetPropType<MomentumTypeTag, Properties::GridGeometry>;
    using MassGridGeometry = GetPropType<MassTypeTag, Properties::GridGeometry>;
    using BasicGridGeometry = typename MassGridGeometry::BasicGridGeometry;

    auto basicGridGeometry  = std::make_shared<BasicGridGeometry>(leafGridView);
    auto momentumGridGeometry = std::make_shared<MomentumGridGeometry>(basicGridGeometry);
    auto massGridGeometry = std::make_shared<MassGridGeometry>(basicGridGeometry);

    using CouplingManager = GetPropType<MomentumTypeTag, Properties::CouplingManager>;
    auto couplingManager = std::make_shared<CouplingManager>(momentumGridGeometry, massGridGeometry);

    using MomentumProblem = GetPropType<MomentumTypeTag, Properties::Problem>;
    using MassProblem = GetPropType<MassTypeTag, Properties::Problem>;
    auto momentumProblem = std::make_shared<MomentumProblem>(momentumGridGeometry, couplingManager);
    auto massProblem = std::make_shared<MassProblem>(massGridGeometry, couplingManager);

    constexpr auto momentumIdx = CouplingManager::momentumIdx;
    constexpr auto massIdx = CouplingManager::massIdx;

    using SolutionVector = typename MDTraits::SolutionVector;
    SolutionVector x;
    x[momentumIdx].resize(momentumGridGeometry->numDofs());
    x[massIdx].resize(massGridGeometry->numDofs());

    using MomentumGridVariables = GetPropType<MomentumTypeTag, Properties::GridVariables>;
    using MassGridVariables = GetPropType<MassTypeTag, Properties::GridVariables>;
    auto momentumGridVariables = std::make_shared<MomentumGridVariables>(momentumProblem, momentumGridGeometry);
    auto massGridVariables = std::make_shared<MassGridVariables>(massProblem, massGridGeometry);

    std::cout << "Initializing coupling manager and grid variables..." << std::endl;

    couplingManager->init(momentumProblem, massProblem, x);
    momentumGridVariables->init(x[momentumIdx]);
    massGridVariables->init(x[massIdx]);

    using Scalar = GetPropType<MomentumTypeTag, Properties::Scalar>;
    const auto tEnd   = getParam<Scalar>("TimeLoop.TEnd");
    const auto dtInit = getParam<Scalar>("TimeLoop.DtInitial");
    const auto dtMax  = getParam<Scalar>("TimeLoop.MaxTimeStepSize", tEnd);
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0.0, dtInit, tEnd);
    timeLoop->setMaxTimeStepSize(dtMax);

    using MomentumSolutionVector = typename MDTraits::template SubDomain<momentumIdx>::SolutionVector;
    auto newmarkBeta = std::make_shared<Experimental::NewmarkBeta<Scalar, MomentumSolutionVector>>();
    newmarkBeta->initialize(x[momentumIdx]);
    momentumProblem->setNewmarkScheme(newmarkBeta);

    VtkOutputModule<MomentumGridVariables, typename MDTraits::template SubDomain<0>::SolutionVector>
        momVtkWriter(*momentumGridVariables, x[momentumIdx], "gripper_momentum");
    momVtkWriter.addVolumeVariable([](const auto& v){
        return Dune::FieldVector<double, 2>{v.displacement(0), v.displacement(1)};
    }, "d");

    std::vector<int> marker(momentumGridGeometry->gridView().size(0));
    for (const auto& element : elements(momentumGridGeometry->gridView()))
        marker[momentumGridGeometry->elementMapper().index(element)] =
            momentumProblem->spatialParams().isActiveLayer(element) ? 1 : 0;

    momVtkWriter.addField(marker, "activeLayer");
    momVtkWriter.write(0.0);

    VtkOutputModule<MassGridVariables, typename MDTraits::template SubDomain<1>::SolutionVector>
        massVtkWriter(*massGridVariables, x[massIdx], "gripper_mass");
    massVtkWriter.addVolumeVariable([](const auto& v){ return v.pressure(); }, "p");
    massVtkWriter.addVolumeVariable([](const auto& v){ return v.solidBulkPressure(); }, "p_solid");
    massVtkWriter.write(0.0);

    auto xOld = x;

    using Assembler = MultiDomainFVAssembler<MDTraits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(
        std::make_tuple(momentumProblem, massProblem),
        std::make_tuple(momentumGridGeometry, massGridGeometry),
        std::make_tuple(momentumGridVariables, massGridVariables),
        couplingManager,
        timeLoop, xOld /* reference previous solution */
    );

    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();

    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    const int vtkInterval = getParam<int>("Vtk.OutputEvery", 1);

    timeLoop->start(); do
    {
        const Scalar dtStep = timeLoop->timeStepSize();

        // Update mass problem with new time (controls swelling BC).
        massProblem->setTime(timeLoop->time() + dtStep);

        nonLinearSolver.solve(x, *timeLoop);

        // Update dynamic time integration state after successful solve.
        newmarkBeta->update(dtStep, x[momentumIdx]);

        // ---- advance ----
        xOld = x;
        momentumGridVariables->advanceTimeStep();
        massGridVariables->advanceTimeStep();

        // Save J^n for the mass storage term in the next time step.
        couplingManager->savePrevTimeStepJacobian();

        timeLoop->advanceTimeStep();

        if (timeLoop->timeStepIndex() % vtkInterval == 0 || timeLoop->finished())
        {
            momVtkWriter.write(timeLoop->time());
            massVtkWriter.write(timeLoop->time());
        }

        timeLoop->reportTimeStep();
        timeLoop->setTimeStepSize(
            nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());
    return 0;
}
