// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief 3D Channel flow test for the staggered grid (Navier-)Stokes model
 */

#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_sub.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include <dumux/freeflow/navierstokes/velocityoutput.hh>
#include <dumux/freeflow/navierstokes/fluxoveraxisalignedsurface.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using MomentumTypeTag = Properties::TTag::ThreeDChannelTestMomentum;
    using MassTypeTag = Properties::TTag::ThreeDChannelTestMass;

    // maybe initialize MPI and/or multithreading backend
    initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // create a grid
    using Scalar = GetPropType<MassTypeTag, Properties::Scalar>;
    using Grid = GetPropType<MassTypeTag, Properties::Grid>;
    Dumux::GridManager<Grid> gridManager;

#if HAVE_DUNE_SUBGRID && GRID_DIM == 3
    const bool isStaircaseGeometry = getParam<bool>("Problem.IsStaircaseGeometry", false);

    auto selector = [&](const auto& element)
    {
        if (!isStaircaseGeometry)
            return true;

        const Scalar deltaX = 0.003;
        const Scalar deltaZ = 0.000075;
        const Scalar deltaY = 0.0003;

        const Scalar eps = 1e-8;
        const auto globalPos = element.geometry().center();

        return globalPos[2] > (deltaZ/deltaX * globalPos[0] + deltaZ/deltaY * globalPos[1] - deltaZ + eps);
    };

    gridManager.init(selector, "Internal");
#elif HAVE_DUNE_SUBGRID
    gridManager.init([&](const auto& element) { return true; });
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

    // the problems (boundary conditions)
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

    // compute coupling stencil and afterwards initialize grid variables (need coupling information)
    couplingManager->init(momentumProblem, massProblem, std::make_tuple(momentumGridVariables, massGridVariables), x);
    massGridVariables->init(x[massIdx]);
    momentumGridVariables->init(x[momentumIdx]);

    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(momentumProblem, massProblem),
                                                 std::make_tuple(momentumGridGeometry, massGridGeometry),
                                                 std::make_tuple(momentumGridVariables, massGridVariables),
                                                 couplingManager);

    // initialize the vtk output module
    using IOFields = GetPropType<MassTypeTag, Properties::IOFields>;
    VtkOutputModule vtkWriter(*massGridVariables, x[massIdx], massProblem->name());
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.addVelocityOutput(std::make_shared<NavierStokesVelocityOutput<MassGridVariables>>());
    vtkWriter.write(0.0);

    // the linear solver
    using LinearSolver = Dumux::UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    // set up three planes over which fluxes are calculated
    FluxOverAxisAlignedSurface flux(*massGridVariables, x[massIdx], assembler->localResidual(massIdx));

    using GridView = typename GetPropType<MassTypeTag, Properties::GridGeometry>::GridView;
    using GlobalPosition = Dune::FieldVector<Scalar, GridView::dimensionworld>;

    const Scalar xMin = massGridGeometry->bBoxMin()[0];
    const Scalar xMax = massGridGeometry->bBoxMax()[0];
    const Scalar yMin = massGridGeometry->bBoxMin()[1];
    const Scalar yMax = massGridGeometry->bBoxMax()[1];

#if GRID_DIM == 3
    const Scalar zMin = massGridGeometry->bBoxMin()[2];
    const Scalar zMax = massGridGeometry->bBoxMax()[2];
#endif

    // the first plane is at the inlet
#if GRID_DIM == 3
    const auto inletLowerLeft = GlobalPosition{xMin, yMin, zMin};
    const auto inletUpperRight = GlobalPosition{xMin, yMax, zMax};
    flux.addAxisAlignedSurface("inlet", inletLowerLeft, inletUpperRight);
#else
    const auto inletLowerLeft = GlobalPosition{xMin, yMin};
    const auto inletUpperRight = GlobalPosition{xMin, yMax};
    flux.addAxisAlignedSurface("inlet", inletLowerLeft, inletUpperRight);
#endif

    // the second plane is at the middle of the channel
    const Scalar planePosMiddleX = xMin + 0.5*(xMax - xMin);
#if GRID_DIM == 3
    const auto middleLowerLeft = GlobalPosition{planePosMiddleX, yMin, zMin};
    const auto middleUpperRight = GlobalPosition{planePosMiddleX, yMax, zMax};
    flux.addAxisAlignedSurface("middle", middleLowerLeft, middleUpperRight);
#else
    const auto middleLowerLeft = GlobalPosition{planePosMiddleX, yMin};
    const auto middleUpperRight = GlobalPosition{planePosMiddleX, yMax};
    flux.addAxisAlignedSurface("middle", middleLowerLeft, middleUpperRight);
#endif

    // The last plane is placed at the outlet of the channel.
#if GRID_DIM == 3
    const auto outletLowerLeft = GlobalPosition{xMax, yMin, zMin};
    const auto outletUpperRight = GlobalPosition{xMax, yMax, zMax};
    flux.addAxisAlignedSurface("outlet", outletLowerLeft, outletUpperRight);
#else
    const auto outletLowerLeft = GlobalPosition{xMax, yMin};
    const auto outletUpperRight = GlobalPosition{xMax, yMax};
    flux.addAxisAlignedSurface("outlet", outletLowerLeft, outletUpperRight);
#endif

    // linearize & solve
    Dune::Timer timer;
    nonLinearSolver.solve(x);

    // write vtk output
    vtkWriter.write(1.0);

    // calculate and print mass fluxes over the planes
    flux.calculateAllFluxes();
    std::cout << "mass flux at inlet is: " << flux.flux("inlet") << std::endl;
    std::cout << "mass flux at middle is: " << flux.flux("middle") << std::endl;
    std::cout << "mass flux at outlet is: " << flux.flux("outlet") << std::endl;
    std::cout << "analyticalFlux: " << massProblem->analyticalFlux()*1e3 << std::endl;

    timer.stop();

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
