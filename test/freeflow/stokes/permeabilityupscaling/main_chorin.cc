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
#include <dune/istl/io.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/io/container.hh>
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

#include "properties_chorin.hh"

#include <dumux/linear/pdesolver.hh>
#include <dumux/io/grid/porenetwork/gridmanager.hh>
#include <dumux/porenetwork/common/pnmvtkoutputmodule.hh>
#include <dumux/porenetwork/common/boundaryflux.hh>
#include <dumux/assembly/fvassembler.hh>

#include "pnm/properties.hh"

std::vector<double> solvePNM()
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::PNMUpscaling;

    // grid
    using GridManager = PoreNetwork::GridManager<3>;
    GridManager gridManager;
    gridManager.init("PNM");
    const auto& leafGridView = gridManager.grid().leafGridView();
    auto gridData = gridManager.getGridData();
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView, *gridData);

    // problem
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;
    auto spatialParams = std::make_shared<SpatialParams>(gridGeometry);
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry, spatialParams, "PNM");

    // instantiate and initialize the discrete and exact solution vectors
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    using VtkOutputFields = GetPropType<TypeTag, Properties::IOFields>;
    using VtkWriter = PoreNetwork::VtkOutputModule<GridVariables, GetPropType<TypeTag, Properties::FluxVariables>, SolutionVector>;
    VtkWriter vtkWriter(*gridVariables, x, problem->name());
    VtkOutputFields::initOutputModule(vtkWriter);
    vtkWriter.addField(gridGeometry->poreVolume(), "poreVolume", Vtk::FieldType::vertex);
    vtkWriter.addField(gridGeometry->throatShapeFactor(), "throatShapeFactor", Vtk::FieldType::element);
    vtkWriter.addField(gridGeometry->throatCrossSectionalArea(), "throatCrossSectionalArea", Vtk::FieldType::element);
    std::vector<int> regionLabel(gridGeometry->numDofs(), -1);
    for (const auto& vertex : vertices(leafGridView))
    {
        const auto vIdx = leafGridView.indexSet().index(vertex);
        const auto& params = gridData->parameters(vertex);
        regionLabel[vIdx] = params[3];
    }
    vtkWriter.addField(regionLabel, "regionLabel", Vtk::FieldType::vertex);

    // assemble & solve
    using Assembler = FVAssembler<TypeTag, DiffMethod::analytic>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);
    using LinearSolver = UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();
    LinearPDESolver<Assembler, LinearSolver> solver(assembler,  linearSolver);
    solver.setVerbose(false);

    using GlobalPosition = typename GridGeometry::GlobalCoordinate;
    const auto sideLengths = getParam<GlobalPosition>("Problem.SideLength");
    problem->setSideLengths(sideLengths); // for setting correct pressure gradient
    problem->setDirection(0); // x-direction
    solver.solve(x);
    vtkWriter.write(0);

    const int maxLabel = *std::max_element(regionLabel.begin(), regionLabel.end());
    std::vector<double> pnmPressures(maxLabel + 1, 0.0);
    for (const auto& vertex : vertices(leafGridView))
    {
        const auto vIdx = leafGridView.indexSet().index(vertex);
        pnmPressures[regionLabel[vIdx]] = x[vIdx][0];
    }

    return pnmPressures;
}

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
    gridManager.grid().globalRefine(getParam<int>("Grid.RefineLevel", 0));

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

    // read regions (map from host grid to pore space grid)
    const std::string regionFileName = getParam<std::string>("Grid.RegionFile");
    const auto regionMarkerHost = readFileToContainer<std::vector<int>>(regionFileName);
    std::vector<int> regionMarker(leafGridView.size(0), -1);
    for (const auto& element : elements(leafGridView))
    {
        const auto eIdx = massGridGeometry->elementMapper().index(element);
        const auto hostEIdx = gridManager.grid().getHostGrid().leafGridView().indexSet().index(
            gridManager.grid().getHostEntity<0>(element)
        );
        regionMarker[eIdx] = regionMarkerHost[hostEIdx];
    }

    // the solution vector
    constexpr auto momentumIdx = CouplingManager::freeFlowMomentumIndex;
    constexpr auto massIdx = CouplingManager::freeFlowMassIndex;
    using SolutionVector = typename Traits::SolutionVector;
    SolutionVector x;
    x[momentumIdx].resize(momentumGridGeometry->numDofs());
    x[massIdx].resize(massGridGeometry->numDofs());
    x = 0.0;

    // initial guess for pressure with PNM
    if (getParam<bool>("Problem.UsePNMInitialPressure", false))
    {
        auto pnmPressures = solvePNM();
        for (const auto& element : elements(leafGridView))
        {
            const auto eIdx = massGridGeometry->elementMapper().index(element);
            x[massIdx][eIdx] = pnmPressures[regionMarker[eIdx]];
        }
    }

    // old solution
    auto xOld = x;

    // the grid variables
    using MomentumGridVariables = GetPropType<MomentumTypeTag, Properties::GridVariables>;
    auto momentumGridVariables = std::make_shared<MomentumGridVariables>(momentumProblem, momentumGridGeometry);
    using MassGridVariables = GetPropType<MassTypeTag, Properties::GridVariables>;
    auto massGridVariables = std::make_shared<MassGridVariables>(massProblem, massGridGeometry);

    // compute coupling stencil and afterwards initialize grid variables (need coupling information)
    couplingManager->init(momentumProblem, massProblem, std::make_tuple(momentumGridVariables, massGridVariables), x, xOld);
    massGridVariables->init(x[massIdx]);
    momentumGridVariables->init(x[momentumIdx]);

    // fake time-dependent problem
    const double dt = getParam<double>("Problem.CFL")*getParam<std::vector<double>>("Grid.PixelDimensions")[0]
        *std::sqrt(getParam<double>("Problem.ArtificialCompressibility"));

    // instantiate time loop
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, dt*1e5);

    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric, /*implicit*/false>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(momentumProblem, massProblem),
                                                 std::make_tuple(momentumGridGeometry, massGridGeometry),
                                                 std::make_tuple(momentumGridVariables, massGridVariables),
                                                 couplingManager, timeLoop, xOld);

    // initialize the vtk output module
    using IOFields = GetPropType<MassTypeTag, Properties::IOFields>;
    VtkOutputModule vtkWriter(*massGridVariables, x[massIdx], massProblem->name());
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.addVelocityOutput(std::make_shared<NavierStokesVelocityOutput<MassGridVariables>>());
    vtkWriter.addField(regionMarker, "region");
    vtkWriter.write(0.0, Dune::VTK::appendedraw);

    using LinearSolver = Dumux::ExplicitDiagonalSolver;
    auto linearSolver = std::make_shared<LinearSolver>();

    LinearPDESolver<Assembler, LinearSolver> solver(assembler,  linearSolver);
    int outputInterval = getParam<int>("Vtk.OutputInterval", 10);

    // assemble and reuse matrix
    assembler->assembleJacobianAndResidual(x);
    solver.reuseMatrix(true);

    if (getParam<bool>("LinearSolver.PrintMatrix", false))
    {
        const auto& jac = assembler->jacobian();
        Dune::printmatrix(std::cout, jac[momentumIdx][momentumIdx], "", "");
        Dune::printmatrix(std::cout, jac[momentumIdx][massIdx], "", "");
        Dune::printmatrix(std::cout, jac[massIdx][momentumIdx], "", "");
        Dune::printmatrix(std::cout, jac[massIdx][massIdx], "", "");

        const auto& res = assembler->residual();
        Dune::printvector(std::cout, res[momentumIdx], "", "");
        Dune::printvector(std::cout, res[massIdx], "", "");
    }

    // time loop
    timeLoop->start(); do
    {
        // solve the linear system
        solver.solve(x);

        // make the new solution the old solution
        xOld = x;
        couplingManager->updateSolution(x);
        momentumGridVariables->advanceTimeStep();
        massGridVariables->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // write vtk output
        if (!(timeLoop->timeStepIndex() % outputInterval))
            vtkWriter.write(timeLoop->time(), Dune::VTK::appendedraw);

        // report statistics of this time step
        timeLoop->reportTimeStep();

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

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
    if (!getParam<bool>("Test.DisableChecks", false))
    {
        const auto KRef = 3.43761e-13;
        const auto phiRef = 0.439627;
        if (Dune::FloatCmp::ne(K, KRef, 1e-5))
            DUNE_THROW(Dune::Exception, "Permeability " << K << " doesn't match reference " << KRef << " by " << K-KRef);
        if (Dune::FloatCmp::ne(phi, phiRef, 1e-5))
            DUNE_THROW(Dune::Exception, "Porosity " << phi << " doesn't match reference " << phiRef << " by " << phi-phiRef);
    }

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
