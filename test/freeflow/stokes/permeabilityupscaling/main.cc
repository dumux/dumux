// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/

#include <config.h>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/common/float_cmp.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_sub.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#include <dumux/linear/incompressiblestokessolver.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/linear/linearsolverparameters.hh>
#include <dumux/linear/linearsolvertraits.hh>

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
    gridManager.init();

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using MomentumGridGeometry = GetPropType<MomentumTypeTag, Properties::GridGeometry>;
    auto momentumGridGeometry = std::make_shared<MomentumGridGeometry>(leafGridView);
    using MassGridGeometry = GetPropType<MassTypeTag, Properties::GridGeometry>;
    auto massGridGeometry = std::make_shared<MassGridGeometry>(leafGridView);

    // the coupling manager
    using Traits = MultiDomainTraits<MomentumTypeTag, MassTypeTag>;
    using CouplingManager = StaggeredFreeFlowCouplingManager<Traits>;
    auto couplingManager = std::make_shared<CouplingManager>();

    // the problems (boundary conditions)
    using MomentumProblem = GetPropType<MomentumTypeTag, Properties::Problem>;
    auto momentumProblem = std::make_shared<MomentumProblem>(momentumGridGeometry, couplingManager);
    using MassProblem = GetPropType<MassTypeTag, Properties::Problem>;
    auto massProblem = std::make_shared<MassProblem>(massGridGeometry, couplingManager);

    // the solution vector
    constexpr auto momentumIdx = CouplingManager::freeFlowMomentumIndex;
    constexpr auto massIdx = CouplingManager::freeFlowMassIndex;
    using SolutionVector = typename Traits::SolutionVector;
    SolutionVector x;
    x[momentumIdx].resize(momentumGridGeometry->numDofs());
    x[massIdx].resize(massGridGeometry->numDofs());
    x = 0.0;

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

    // solve
    if (getParam<bool>("LinearSolver.UseDirectSolver", false))
    {
        using LinearSolver = UMFPackBackend;
        auto linearSolver = std::make_shared<LinearSolver>();
        // the non-linear solver
        using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
        NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);
        nonLinearSolver.solve(x);
    }
    else if (getParam<bool>("LinearSolver.UseUzawa", false))
    {
        using LinearSolver = UzawaBiCGSTABBackend<LinearSolverTraits<MassGridGeometry>>;
        auto linearSolver = std::make_shared<LinearSolver>();
        // the non-linear solver
        using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
        NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);
        nonLinearSolver.solve(x);
    }
    else
    {
        using LinearSolver = IncompressibleStokesSolver<typename Assembler::JacobianMatrix, typename Assembler::ResidualType>;
        auto params = LinearSolverParameters<LinearSolverTraits<MassGridGeometry>>::createParameterTree();
        params["Component.LiquidDensity"] = getParam<std::string>("Component.LiquidDensity");
        params["Component.LiquidKinematicViscosity"] = getParam<std::string>("Component.LiquidKinematicViscosity");
        params["preconditioner.UseBlockTriangular"] = getParam<std::string>("LinearSolver.Preconditioner.UseBlockTriangular", "false");
        params["preconditioner.UseFullSchurComplement"] = getParam<std::string>("LinearSolver.Preconditioner.UseFullSchurComplement", "false");
        params["preconditioner.UseDirectVelocitySolver"] = getParam<std::string>("LinearSolver.Preconditioner.UseDirectVelocitySolver", "false");
        const auto [massMatrix, dirichletDofs]
            = computeIncompressibleStokesMassMatrixAndDirichletDofs<SolutionVector>(
                leafGridView,
                *massGridGeometry, massIdx,
                *momentumGridGeometry, *momentumProblem, momentumIdx
            );
        auto linearSolver = std::make_shared<LinearSolver>(massMatrix, dirichletDofs, params);
        using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
        NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);
        nonLinearSolver.solve(x);
    }

    // post-processing and output
    vtkWriter.write(1.0);

    // set up two planes over which fluxes are calculated
    FluxOverAxisAlignedSurface flux(*massGridVariables, x[massIdx], assembler->localResidual(massIdx));

    using GridView = typename GetPropType<MassTypeTag, Properties::GridGeometry>::GridView;
    using GlobalPosition = Dune::FieldVector<Scalar, GridView::dimensionworld>;

    const Scalar xMin = massGridGeometry->bBoxMin()[0];
    const Scalar xMax = massGridGeometry->bBoxMax()[0];
    const Scalar yMin = massGridGeometry->bBoxMin()[1];
    const Scalar yMax = massGridGeometry->bBoxMax()[1];
    const Scalar zMin = massGridGeometry->bBoxMin()[2];
    const Scalar zMax = massGridGeometry->bBoxMax()[2];

    // the first plane is at the inlet
    const auto inletLowerLeft = GlobalPosition{xMin, yMin, zMin};
    const auto inletUpperRight = GlobalPosition{xMin, yMax, zMax};
    flux.addAxisAlignedSurface("inlet", inletLowerLeft, inletUpperRight);

    // The last plane is placed at the outlet of the channel.
    const auto outletLowerLeft = GlobalPosition{xMax, yMin, zMin};
    const auto outletUpperRight = GlobalPosition{xMax, yMax, zMax};
    flux.addAxisAlignedSurface("outlet", outletLowerLeft, outletUpperRight);

    // calculate and print mass fluxes over the planes
    flux.calculateAllFluxes();
    std::cout << "mass flux at inlet is: " << flux.flux("inlet") << " µm^3/s" <<  std::endl;
    std::cout << "mass flux at outlet is: " << flux.flux("outlet") << " µm^3/s" << std::endl;

    // make sure solver is converged enough
    if (Dune::FloatCmp::ne(flux.flux("inlet"), -flux.flux("outlet"), 1e-5))
        DUNE_THROW(Dune::Exception, "Inlet and outlet fluxes do not match by " << flux.flux("inlet")+flux.flux("outlet"));

    const auto cells = getParam<std::array<double, Grid::dimension>>("Grid.Cells");
    const auto dims = getParam<std::array<double, Grid::dimension>>("Grid.PixelDimensions");
    const auto totalAreaOutlet = dims[1]*cells[1]*dims[2]*cells[2];
    // viscosity, density, and pressure gradient are 1.0
    // With the pressure difference across the domain defined as 1/m,
    // the permeability can be determined by dividing the flux at the outlet by the area of the outlet.
    const auto K = std::abs(flux.flux("outlet")/totalAreaOutlet)*1e-12; // domain is in µm
    // The porosity can be calculated for uniform sized grid cells by dividing
    // the number of cells in the leafGridView by the number of cells used in the original hostgrid.
    const auto phi = double(leafGridView.size(0))/double(cells[0]*cells[1]*cells[2]);
    std::cout << "Permeability is: " << K << " m^2" << std::endl;
    std::cout << "Porosity is: " << phi << std::endl;

    // test against reference solution (hard-coded here)
    // note that for the test the geometry is very coarse to minimize runtime
    // therefore these reference values are just regression test references
    // computed with this test at the time it was set up
    const auto KRef = 3.43761e-13;
    const auto phiRef = 0.439627;
    if (Dune::FloatCmp::ne(K, KRef, 1e-5))
        DUNE_THROW(Dune::Exception, "Permeability " << K << " doesn't match reference " << KRef << " by " << K-KRef);
    if (Dune::FloatCmp::ne(phi, phiRef, 1e-5))
        DUNE_THROW(Dune::Exception, "Porosity " << phi << " doesn't match reference " << phiRef << " by " << phi-phiRef);

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
