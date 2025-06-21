// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Stationary test for the staggered grid Navier-Stokes model with periodic BC
 */

#ifndef MOMENTUM_TYPETAG
#define MOMENTUM_TYPETAG PeriodicTestMomentumStaggered
#endif

#ifndef MASS_TYPETAG
#define MASS_TYPETAG PeriodicTestMassTpfa
#endif

#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/io/grid/gridmanager_sp.hh>
#include <dumux/io/grid/gridmanager_sub.hh>
#include <dumux/io/grid/periodicgridtraits.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/vtk/function.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>

#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include <dumux/io/vtk/intersectionwriter.hh>
#include <dumux/freeflow/navierstokes/velocityoutput.hh>

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

        if constexpr (GridGeometry::discMethod == Dumux::DiscretizationMethods::fcstaggered)
        {
            for (const auto& scv : scvs(fvGeometry))
            {
                const auto& vars = elemVolVars[scv];
                faceVelocity[scv.dofIndex()][scv.dofAxis()] = vars.velocity();
            }
        }
        else if constexpr (GridGeometry::discMethod == Dumux::DiscretizationMethods::fcdiamond)
        {
            const auto elemGeo = element.geometry();
            const auto elemSol = elementSolution(element, x, gridGeometry);
            for (const auto& scv : scvs(fvGeometry))
                faceVelocity[scv.dofIndex()] = elemVolVars[scv].velocity();
        }
        else
            DUNE_THROW(Dune::Exception, "Invalid discretization type: " << GridGeometry::discMethod);
    }
}

int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using MomentumTypeTag = Properties::TTag::MOMENTUM_TYPETAG;
    using MassTypeTag = Properties::TTag::MASS_TYPETAG;

    // maybe initialize MPI and/or multithreading backend
    initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // try to create a grid (from the given grid file or the input file)
    using Grid = GetPropType<MomentumTypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;

#if HAVE_DUNE_SUBGRID && USESUBGRID
    auto selector = [](const auto& element) {
        const static auto trivialSelector = getParam<bool>("Problem.TrivialSelector", true);
        const static auto radius = getParam<double>("Problem.Radius", 0.25);
        const auto distance = std::max(std::abs(element.geometry().center()[0] - 0.5),
                                       std::abs(element.geometry().center()[1] - 0.5));
        return trivialSelector || distance > radius - 1e-6;
    };
    gridManager.init(selector);

    // Verify that the subgrid still has a conforming periodic boundary
    // Take additional care with boundary conditions if not guaranteed
    PeriodicGridTraits<Grid> periodicGridTraits(gridManager.grid());
    periodicGridTraits.verifyConformingPeriodicBoundary();
#else
    gridManager.init();
#endif

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using MomentumGridGeometry = GetPropType<MomentumTypeTag, Properties::GridGeometry>;
    auto momentumGridGeometry = std::make_shared<MomentumGridGeometry>(leafGridView);
    using MassGridGeometry = GetPropType<MassTypeTag, Properties::GridGeometry>;
    auto massGridGeometry = std::make_shared<MassGridGeometry>(leafGridView);

    // the coupling manager
    using CouplingManager = GetPropType<MomentumTypeTag, Properties::CouplingManager>;
    auto couplingManager = std::make_shared<CouplingManager>();

    // the problem (boundary conditions)
    using MomentumProblem = GetPropType<MomentumTypeTag, Properties::Problem>;
    auto momentumProblem = std::make_shared<MomentumProblem>(momentumGridGeometry, couplingManager);
    using MassProblem = GetPropType<MassTypeTag, Properties::Problem>;
    auto massProblem = std::make_shared<MassProblem>(massGridGeometry, couplingManager);

    // the solution vector
    constexpr auto momentumIdx = CouplingManager::freeFlowMomentumIndex;
    constexpr auto massIdx = CouplingManager::freeFlowMassIndex;
    using Traits = MultiDomainTraits<MomentumTypeTag, MassTypeTag>;
    using SolutionVector = typename Traits::SolutionVector;
    SolutionVector x;
    x[momentumIdx].resize(momentumGridGeometry->numDofs());
    x[massIdx].resize(massGridGeometry->numDofs());

    // the grid variables
    using MomentumGridVariables = GetPropType<MomentumTypeTag, Properties::GridVariables>;
    auto momentumGridVariables = std::make_shared<MomentumGridVariables>(momentumProblem, momentumGridGeometry);
    using MassGridVariables = GetPropType<MassTypeTag, Properties::GridVariables>;
    auto massGridVariables = std::make_shared<MassGridVariables>(massProblem, massGridGeometry);

    // initialize coupling stencil first and then grid variables (need coupling variables)
    couplingManager->init(momentumProblem, massProblem, std::make_tuple(momentumGridVariables, massGridVariables), x);
    massGridVariables->init(x[massIdx]);
    momentumGridVariables->init(x[momentumIdx]);

    // the assembler with time loop for instationary problem
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(momentumProblem, massProblem),
                                                 std::make_tuple(momentumGridGeometry, massGridGeometry),
                                                 std::make_tuple(momentumGridVariables, massGridVariables),
                                                 couplingManager);

    // initialize the vtk output module
    using IOFields = GetPropType<MassTypeTag, Properties::IOFields>;
    constexpr bool isDiamond = MassGridGeometry::discMethod == DiscretizationMethods::fcdiamond;
    constexpr auto mode = isDiamond ? Dune::VTK::nonconforming : Dune::VTK::conforming;
    VtkOutputModule vtkWriter(*massGridVariables, x[massIdx], massProblem->name(), "", mode);
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.addVelocityOutput(std::make_shared<NavierStokesVelocityOutput<MassGridVariables>>());
    ConformingIntersectionWriter faceVtk(momentumGridGeometry->gridView());

    using LinearSolver = Dumux::UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    // linearize & solve
    Dune::Timer timer;
    nonLinearSolver.solve(x);

    if constexpr (MomentumGridGeometry::discMethod != DiscretizationMethods::pq1bubble)
    {
        std::vector<Dune::FieldVector<double,2>> faceVelocityVector(x[momentumIdx].size());
        updateVelocities(faceVelocityVector, *momentumGridGeometry, *momentumGridVariables, x[momentumIdx]);

        std::vector<std::size_t> dofIdx(x[momentumIdx].size());
        for (const auto& facet : facets(momentumGridGeometry->gridView()))
        {
            const auto idx = momentumGridGeometry->gridView().indexSet().index(facet);
            dofIdx[idx] = idx;
        }
        faceVtk.addField(dofIdx, "dofIdx");

        faceVtk.addField(faceVelocityVector, "velocityVector");
        faceVtk.write("facedata", Dune::VTK::ascii);
    }

    // write vtk output
    vtkWriter.write(1.0);
    timer.stop();

    // verify periodic mapping of dofs
    {
        const bool periodicInX = getParam<bool>("Problem.PeriodicInX", false);
        const auto eps = 1e-6;
        const auto bBoxMax = momentumGridGeometry->bBoxMax();
        const auto bBoxMin = momentumGridGeometry->bBoxMin();
        auto fvGeometry = localView(*momentumGridGeometry);
        for (const auto& element : elements(momentumGridGeometry->gridView()))
        {
            fvGeometry.bindElement(element);
            for (const auto& scv : scvs(fvGeometry))
            {
                const bool isPeriodic = std::min(scv.dofPosition()[1]-bBoxMin[1],
                                                 bBoxMax[1]-scv.dofPosition()[1]) < eps
                    || (periodicInX && std::min(scv.dofPosition()[0]-bBoxMin[0],
                                                 bBoxMax[0]-scv.dofPosition()[0]) < eps);
                if (isPeriodic != momentumGridGeometry->dofOnPeriodicBoundary(scv.dofIndex()))
                    DUNE_THROW(Dune::Exception, "Grid does not exhibit expected periodicity");
                if (isPeriodic)
                {
                    //const auto periodicScv = fvGeometry.outsidePeriodicScv(scv);
                    const auto periodicallyMappedDofs =
                        Dumux::Deprecated::rangeOfPeriodicallyMappedDofs(*momentumGridGeometry, scv.dofIndex());
                    for (const auto periodicDof : periodicallyMappedDofs)
                    {
                        //const auto distance = scv.dofPosition() - periodicScv.dofPosition();
                        //if ((std::abs(distance[0]) > eps
                        //        || std::abs(std::abs(distance[1]) - (bBoxMax[1]-bBoxMin[1]) ) > eps )
                        //   && (!periodicInX || std::abs(distance[1]) > eps
                        //        || std::abs(std::abs(distance[0]) - (bBoxMax[0]-bBoxMin[0]) ) > eps ))
                        //    DUNE_THROW(Dune::Exception, "Grid does not exhibit expected periodicity");
                        for (int i = 0; i < x[momentumIdx][periodicDof].size(); ++i)
                            if (std::abs(x[momentumIdx][periodicDof][i] - x[momentumIdx][scv.dofIndex()][i]) > eps)
                                DUNE_THROW(Dune::Exception, "Periodicity constraints not enforced");
                    }
                }

            }
        }
    }

    const auto& comm = Dune::MPIHelper::getCommunication();
    std::cout << "Simulation took " << timer.elapsed() << " seconds on "
              << comm.size() << " processes.\n"
              << "The cumulative CPU time was " << timer.elapsed()*comm.size() << " seconds.\n";

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
}
