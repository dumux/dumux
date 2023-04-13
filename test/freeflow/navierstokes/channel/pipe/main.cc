// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/float_cmp.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/partial.hh>

#include <dumux/geometry/diameter.hh>

#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/vtkoutputmodule.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>

#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include <dumux/freeflow/navierstokes/velocityoutput.hh>
#include <dumux/freeflow/navierstokes/fluxoveraxisalignedsurface.hh>
#include <test/freeflow/navierstokes/analyticalsolutionvectors.hh>
#include <test/freeflow/navierstokes/errors.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using MomentumTypeTag = Properties::TTag::PipeFlowMomentum;
    using MassTypeTag = Properties::TTag::PipeFlowMass;

    // maybe initialize MPI and/or multithreading backend
    initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // try to create a grid (from the given grid file or the input file)
    Dumux::GridManager<GetPropType<MomentumTypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    // we compute on the leaf grid view
    const auto& gridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using MomentumGridGeometry = GetPropType<MomentumTypeTag, Properties::GridGeometry>;
    auto momentumGridGeometry = std::make_shared<MomentumGridGeometry>(gridView);

    using MassGridGeometry = GetPropType<MassTypeTag, Properties::GridGeometry>;
    auto massGridGeometry = std::make_shared<MassGridGeometry>(gridView);

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
    momentumProblem->applyInitialSolution(x[momentumIdx]);
    massProblem->applyInitialSolution(x[massIdx]);


    // the grid variables
    using MomentumGridVariables = GetPropType<MomentumTypeTag, Properties::GridVariables>;
    auto momentumGridVariables = std::make_shared<MomentumGridVariables>(momentumProblem, momentumGridGeometry);

    using MassGridVariables = GetPropType<MassTypeTag, Properties::GridVariables>;
    auto massGridVariables = std::make_shared<MassGridVariables>(massProblem, massGridGeometry);

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

    nonLinearSolver.solve(x);
    vtkWriter.write(1.0);

    // the discrete L2 and Linfity errors
    const bool printErrors = getParam<bool>("Problem.PrintErrors", false);
    const bool printConvergenceTestFile = getParam<bool>("Problem.PrintConvergenceTestFile", false);
    NavierStokesTest::Errors errors(momentumProblem, massProblem, x);
    NavierStokesTest::ErrorCSVWriter errorCSVWriter(momentumProblem, massProblem);

    if (printErrors || printConvergenceTestFile)
    {
        NavierStokesTest::Errors errors(momentumProblem, massProblem, x);
        NavierStokesTest::ErrorCSVWriter errorCSVWriter(momentumProblem, massProblem);
        errorCSVWriter.printErrors(errors);

        if (printConvergenceTestFile)
            convergenceTestAppendErrors(momentumProblem, massProblem, errors);
    }

    // Use FluxOverAxisAlignedSurface class to calculate fluxes entering and leaving the domain.
    auto fluxOverSurface = FluxOverAxisAlignedSurface(*massGridVariables, x[massIdx], assembler->localResidual(massIdx));
    using GlobalPosition = typename MassGridGeometry::SubControlVolume::GlobalPosition;

    // Add a surface at the lower bottom (inlet).
    GlobalPosition inletLowerLeft = massGridGeometry->bBoxMin();
    GlobalPosition inletUpperRight = GlobalPosition{
        massGridGeometry->bBoxMax()[0], massGridGeometry->bBoxMin()[1]
    };

    // The surface has to coincide with the element's faces. Let's assume the user does not give that position
    // correctly but places the surface somewhere in-between. The FluxOverAxisAlignedSurface class
    // will automatically shift that surface to the closest valid position.
    const auto firstElement = (*elements(massGridGeometry->gridView()).begin());
    const auto offset = 0.25*diameter(firstElement.geometry());
    inletLowerLeft[1] += offset;
    inletUpperRight[1] += offset;
    fluxOverSurface.addAxisAlignedSurface("inlet", inletLowerLeft, inletUpperRight);

    // Add a plane (infinite extend) at the outlet.
    const GlobalPosition planeCenter = massGridGeometry->bBoxMax(); // the exact position in x-direction does not matter
    const std::size_t directionIndex = 1; // plane normal points in y-direction
    fluxOverSurface.addAxisAlignedPlane("outlet", planeCenter, directionIndex);

    // Calculate and print fluxes
    fluxOverSurface.calculateAllFluxes();
    fluxOverSurface.printAllFluxes();

    // Sanity check
    if (Dune::FloatCmp::ne(fluxOverSurface.flux("inlet"), -fluxOverSurface.flux("outlet")))
        DUNE_THROW(Dune::InvalidStateException, "Inlet and outlet fluxes differ");

    const auto analyticalMassFlux = massProblem->analyticalMassFlux();
    std::cout << "Analytical mass flux " << analyticalMassFlux << std::endl;

    if (Dune::FloatCmp::ne(fluxOverSurface.flux("outlet")[0], analyticalMassFlux, 5e-3))
        DUNE_THROW(Dune::InvalidStateException, "Outlet flux differs too much from analytical solution");

    // print dumux end message
    if (mpiHelper.rank() == 0)
        Parameters::print();

    return 0;

} // end main
