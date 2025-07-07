// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

 // ## The main file (`main.cc`)
 // [[content]]
 //
 // ### Included header files
 // [[details]] includes
 // [[exclude]]
 // Some generic includes.
#include <config.h>
#include <ctime>
#include <iostream>
// [[/exclude]]

// These are DUNE helper classes related to parallel computations and file I/O.
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>

// The following headers include functionality related to property definition or retrieval, as well as
// the retrieval of input parameters specified in the input file or via the command line.
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/initialize.hh>

// The following files contain the non-linear Newton solver for multi-domain problems, the available linear solver backends and the assembler for the linear
// systems arising from the staggered-grid discretization.
#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/assembly/diffmethod.hh> // analytic or numeric differentiation

// The following class provides a convenient way of writing of dumux simulation results to VTK format and `velocityoutput.hh` allows to additionally write out velocity data of the staggered grid.
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/freeflow/navierstokes/velocityoutput.hh>
// The gridmanager constructs a grid from the information in the input or grid file.
// Many different Dune grid implementations are supported, of which a list can be found
// in `gridmanager.hh`.
#include <dumux/io/grid/gridmanager.hh>

// This class contains functionality for additional flux output.
#include <dumux/freeflow/navierstokes/fluxoveraxisalignedsurface.hh>


// In this header three `TypeTag`s are defined, which collect
// the properties that are required for the simulation.
// It also contains the actual problem with initial and boundary conditions.
// For detailed information, please have a look
// at the documentation provided therein.
#include "properties.hh"
// [[/details]]
//

// ### The main function
// We will now discuss the main program flow implemented within the `main` function.
// At the beginning of each program using Dune, an instance of `Dune::MPIHelper` has to
// be created. Moreover, we parse the run-time arguments from the command line and the
// input file:
// [[codeblock]]
int main(int argc, char** argv) try
{
    using namespace Dumux;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // parse command line arguments and input file
    Parameters::init(argc, argv);
    // [[/codeblock]]

    // We define convenience aliases for the type tags of the subproblems. The type
    // tags contain all the properties that are needed to define the model and the problem
    // setup. Throughout the main file, we will obtain types defined for these type tags
    // using the property system, i.e. with `GetPropType`. Shared properties can be obtained through
    // either of them.
    using MomentumTypeTag = Properties::TTag::ChannelExampleMomentum;
    using MassTypeTag = Properties::TTag::ChannelExampleMass;

    // #### Step 1: Create the grid
    // The `GridManager` class creates the grid from information given in the input file.
    // This can either be a grid file, or in the case of structured grids, one can specify the coordinates
    // of the corners of the grid and the number of cells to be used to discretize each spatial direction.
    // 'Grid` is a property that is shared by both type tags (see `properties.hh`).
    // As stated above, it can therefore be obtained through either `MassTypeTag` of `MomentumTypeTag`.
    // [[codeblock]]
    GridManager<GetPropType<MomentumTypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    // We compute on the leaf grid view.
    const auto& leafGridView = gridManager.grid().leafGridView();
    // [[/codeblock]]

    // #### Step 2: Setting up and solving the problem
    // First, a finite volume grid geometry is constructed from the grid that was created above.
    // This builds the sub-control volumes (`scv`) and sub-control volume faces (`scvf`) for each element
    // of the grid partition.
    // This is done separately for the momentum and mass grid geometries.
    using MomentumGridGeometry = GetPropType<MomentumTypeTag, Properties::GridGeometry>;
    auto momentumGridGeometry = std::make_shared<MomentumGridGeometry>(leafGridView);
    using MassGridGeometry = GetPropType<MassTypeTag, Properties::GridGeometry>;
    auto massGridGeometry = std::make_shared<MassGridGeometry>(leafGridView);

    // We introduce the multidomain coupling manager, which will couple the two subproblems for mass
    // and momentum. The type can be obtained using either of the two type tags.
    using CouplingManager = GetPropType<MomentumTypeTag, Properties::CouplingManager>;
    auto couplingManager = std::make_shared<CouplingManager>();

    // We now instantiate the problems, in which we define the boundary and initial conditions.
    using MassProblem = GetPropType<MassTypeTag, Properties::Problem>;
    auto massProblem = std::make_shared<MassProblem>(massGridGeometry, couplingManager);
    using MomentumProblem = GetPropType<MomentumTypeTag, Properties::Problem>;
    auto momentumProblem = std::make_shared<MomentumProblem>(momentumGridGeometry, couplingManager);

    // We set a solution vector `x` which consist of two parts: one part (indexed by `massIdx`)
    // is for the pressure degrees of freedom (`dofs`) living in grid cell centers. Another part
    // (indexed by `momentumIdx`) is for degrees of freedom defining the normal velocities on grid cell faces.
    // The relevant types can be accessed through the `MultiDomainTraits` of the coupled problem.
    // We initialize the solution vector by what was defined as the initial solution of the problem.
    using Traits = MultiDomainTraits<MomentumTypeTag, MassTypeTag>;
    using SolutionVector = typename Traits::SolutionVector;
    constexpr auto momentumIdx = CouplingManager::freeFlowMomentumIndex;
    constexpr auto massIdx = CouplingManager::freeFlowMassIndex;
    SolutionVector x;
    momentumProblem->applyInitialSolution(x[momentumIdx]);
    massProblem->applyInitialSolution(x[massIdx]);

    // The grid variables are used to store variables (primary and secondary variables) of the two subproblems.
    using MomentumGridVariables = GetPropType<MomentumTypeTag, Properties::GridVariables>;
    auto momentumGridVariables = std::make_shared<MomentumGridVariables>(momentumProblem, momentumGridGeometry);
    using MassGridVariables = GetPropType<MassTypeTag, Properties::GridVariables>;
    auto massGridVariables = std::make_shared<MassGridVariables>(massProblem, massGridGeometry);

    // After initializing the coupling manager the coupling context is set up and the grid variables
    // of the subproblems can be initialized.
    couplingManager->init(momentumProblem, massProblem, std::make_tuple(momentumGridVariables, massGridVariables), x);
    momentumGridVariables->init(x[momentumIdx]);
    massGridVariables->init(x[massIdx]);

    // We then initialize the predefined model-specific VTK output.
    using IOFields = GetPropType<MassTypeTag, Properties::IOFields>;
    VtkOutputModule vtkWriter(*massGridVariables, x[massIdx], massProblem->name());
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.addVelocityOutput(std::make_shared<NavierStokesVelocityOutput<MassGridVariables>>());
    vtkWriter.write(0.0);

    // We create and initialize the `assembler` for the stationary problem.
    // This is where the Jacobian matrix for the Newton solver is assembled.
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(momentumProblem, massProblem),
                                                 std::make_tuple(momentumGridGeometry, massGridGeometry),
                                                 std::make_tuple(momentumGridVariables, massGridVariables),
                                                 couplingManager);

    // We use UMFPack as direct linear solver within each Newton iteration.
    using LinearSolver = Dumux::UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();

    // This example considers a linear problem (incompressible Stokes flow), therefore
    // the non-linear Newton solver is not really necessary.
    // For sake of generality, we nevertheless use it here such that the example can be easily
    // changed to a non-linear problem by switching on the inertia terms in the input file or by choosing a compressible fluid.
    // In the following piece of code we instantiate the non-linear newton solver and let it solve
    // the problem.
    // [[codeblock]]
    // alias for and instantiation of the newton solver
    using NewtonSolver = Dumux::MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);
    // [[/codeblock]]

    // [[details]] calculation of surface fluxes
    //
    // We set up two surfaces over which fluxes are calculated.
    // We determine the extent $`[xMin,xMax] \times [yMin,yMax]`$ of the physical domain.
    // The first surface (added by the first call of `addSurface`) shall be placed at the middle of the channel.
    // The second surface (second call of `addSurface`) is placed at the outlet of the channel.
    // [[codeblock]]
    FluxOverAxisAlignedSurface flux(*massGridVariables, x[massIdx], assembler->localResidual(massIdx));

    using Scalar = typename Traits::Scalar;

    const Scalar xMin = massGridGeometry->bBoxMin()[0];
    const Scalar xMax = massGridGeometry->bBoxMax()[0];
    const Scalar yMin = massGridGeometry->bBoxMin()[1];
    const Scalar yMax = massGridGeometry->bBoxMax()[1];

    const Scalar planePosMiddleX = xMin + 0.5*(xMax - xMin);

    using GridView = typename MassGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    const auto p0middle = GlobalPosition{planePosMiddleX, yMin};
    const auto p1middle = GlobalPosition{planePosMiddleX, yMax};
    flux.addAxisAlignedSurface("middle", p0middle, p1middle);

    const auto p0outlet = GlobalPosition{xMax, yMin};
    const auto p1outlet = GlobalPosition{xMax, yMax};
    flux.addAxisAlignedSurface("outlet", p0outlet, p1outlet);

    using FluxVariables = GetPropType<MassTypeTag, Properties::FluxVariables>;
    auto volumeFlux = [&](const auto& element,
                         const auto& fvGeometry,
                         const auto& elemVolVars,
                         const auto& scvf,
                         const auto& elemFluxVarsCache)
    {
        if (scvf.boundary() && massProblem->boundaryTypes(element, scvf).hasNeumann())
            return massProblem->neumann(element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf)[0]/elemVolVars[scvf.insideScvIdx()].density();
        else
        {
            FluxVariables fluxVars;
            fluxVars.init(*massProblem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
            return fluxVars.getAdvectiveFlux([](const auto& volVars) { return 1.0; });
        }
    };
    // [[/codeblock]]
    // [[/details]]
    //

    // Solve the (potentially non-linear) system.
    // [[codeblock]]
    nonLinearSolver.solve(x);
    // [[/codeblock]]
    // In the following we calculate and print mass and volume fluxes over the planes specified above
    // (you have to click to unfold the code showing how to set up the surface fluxes above).
    // [[codeblock]]
    flux.calculateAllFluxes();
    if (GetPropType<MassTypeTag, Properties::ModelTraits>::enableEnergyBalance())
    {
        std::cout << "mass / energy flux at middle is: " << flux.flux("middle") << std::endl;
        std::cout << "mass / energy flux at outlet is: " << flux.flux("outlet") << std::endl;
    }
    else
    {
        std::cout << "mass flux at middle is: " << flux.flux("middle") << std::endl;
        std::cout << "mass flux at outlet is: " << flux.flux("outlet") << std::endl;
    }

    flux.calculateFluxes(volumeFlux);
    std::cout << "volume flux at middle is: " << flux.flux("middle")[0] << std::endl;
    std::cout << "volume flux at outlet is: " << flux.flux("outlet")[0] << std::endl;
    // [[/codeblock]]

    // #### Final Output
    // We write the VTK output and conclude by printing the dumux end message.
    // [[codeblock]]
    vtkWriter.write(1.0);

    if (mpiHelper.rank() == 0)
        Parameters::print();

    return 0;
} // end main
// [[/codeblock]]
// #### Exception handling
// In this part of the main file we catch and print possible exceptions that could
// occur during the simulation.
// [[details]] exception handler
// [[codeblock]]
// errors related to run-time parameters
catch (Dumux::ParameterException &e)
{
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
}
// errors related to the parsing of Dune grid files
catch (Dune::DGFException & e)
{
    std::cerr << "DGF exception thrown (" << e <<
                 "). Most likely, the DGF file name is wrong "
                 "or the DGF file is corrupted, "
                 "e.g. missing hash at end of file or wrong number (dimensions) of entries."
                 << " ---> Abort!" << std::endl;
    return 2;
}
// generic error handling with Dune::Exception
catch (Dune::Exception &e)
{
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
}
// other exceptions
catch (...)
{
    std::cerr << "Unknown exception thrown! ---> Abort!" << std::endl;
    return 4;
}
// [[/codeblock]]
// [[/details]]
// [[/content]]
