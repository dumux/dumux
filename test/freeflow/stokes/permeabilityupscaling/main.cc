// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Pore flow test for permeability upscaling
 */

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

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/stokes_solver.hh>

#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include <dumux/freeflow/navierstokes/velocityoutput.hh>
#include <dumux/freeflow/navierstokes/fluxoveraxisalignedsurface.hh>

#include <dumux/geometry/intersectionentityset.hh>

#include "properties.hh"

template<class Position>
std::vector<Position> corners(const Position& c, const Position& l)
{
    const Position lowerLeft = c - 0.5*l;
    return {
        lowerLeft                              , lowerLeft + Position({l[0], 0.0, 0.0}),
        lowerLeft + Position({0.0, l[1], 0.0}) , lowerLeft + Position({l[0], l[1], 0.0}),
        lowerLeft + Position({0.0, 0.0, l[2]}) , lowerLeft + Position({l[0], 0.0, l[2]}),
        lowerLeft + Position({0.0, l[1], l[2]}), lowerLeft + Position({l[0], l[1], l[2]})
    };
}


// helper function to generate an REV geometry
template<class ctype, int dim>
auto makeAveragingVolume(const Dune::FieldVector<ctype, dim>& center,
                         const Dune::FieldVector<ctype, dim>& l)
{
    return Dune::MultiLinearGeometry<ctype, dim, dim>{
        Dune::GeometryTypes::cube(dim), corners(center, l)
    };
}

template<class Vector, class MomGG, class MassGG, class MomP, class MomIdx, class MassIdx>
auto dirichletDofs(std::shared_ptr<MomGG> momentumGridGeometry,
                   std::shared_ptr<MassGG> massGridGeometry,
                   std::shared_ptr<MomP> momentumProblem,
                   MomIdx momentumIdx, MassIdx massIdx)
{
    Vector dirichletDofs;
    dirichletDofs[momentumIdx].resize(momentumGridGeometry->numDofs());
    dirichletDofs[massIdx].resize(massGridGeometry->numDofs());
    dirichletDofs = 0.0;

    auto fvGeometry = localView(*momentumGridGeometry);
    for (const auto& element : elements(momentumGridGeometry->gridView()))
    {
        fvGeometry.bind(element);
        for (const auto& scvf : scvfs(fvGeometry))
        {
            if (!scvf.boundary() || !scvf.isFrontal())
                continue;
            const auto& scv = fvGeometry.scv(scvf.insideScvIdx());
            if (scv.boundary())
            {
                const auto bcTypes = momentumProblem->boundaryTypes(element, scvf);
                for (int i = 0; i < bcTypes.size(); ++i)
                    if (bcTypes.isDirichlet(i))
                        dirichletDofs[momentumIdx][scv.dofIndex()][i] = 1.0;
            }
        }
    }

    return dirichletDofs;
}

int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using MomentumTypeTag = Properties::TTag::PoreFlowTestMomentum;
    using MassTypeTag = Properties::TTag::PoreFlowTestMass;

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

    // mass-momentum coupling manager
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

    // the linearize and solve
    if (getParam<bool>("LinearSolver.UseIterativeSolver", true))
    {
        using Matrix = typename Assembler::JacobianMatrix;
        using Vector = typename Assembler::ResidualType;
        using LinearSolver = StokesSolver<Matrix, Vector, MomentumGridGeometry, MassGridGeometry>;
        auto dDofs = dirichletDofs<Vector>(momentumGridGeometry, massGridGeometry, momentumProblem, momentumIdx, massIdx);
        auto linearSolver = std::make_shared<LinearSolver>(momentumGridGeometry, massGridGeometry, dDofs);
        using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
        NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);
        nonLinearSolver.solve(x);
    }
    else
    {
        using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
        auto linearSolver = std::make_shared<LinearSolver>();
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
    // do the same calculations by using volume averaging
    ////////////////////////////////////////////////////////////

    // Calculate element velocities
    std::vector<GlobalPosition> velocity(leafGridView.size(0));
    auto fvGeometry = localView(*massGridGeometry);
    for (const auto& element : elements(leafGridView))
    {
        const auto eIdx = massGridGeometry->elementMapper().index(element);
        fvGeometry.bind(element);
        const auto getFaceVelocity = [&](const auto&, const auto& scvf)
        { return massProblem->faceVelocity(element, fvGeometry, scvf); };

        velocity[eIdx] = StaggeredVelocityReconstruction::cellCenterVelocity(getFaceVelocity, fvGeometry);
    }

    const auto size = massGridGeometry->bBoxMax() - massGridGeometry->bBoxMin();
    const auto center = 0.5*(massGridGeometry->bBoxMax() + massGridGeometry->bBoxMin());
    const auto geo = makeAveragingVolume(center, size);

    using GES = GridViewGeometricEntitySet<GridView>;
    using AES = SingleGeometryEntitySet<std::decay_t<decltype(geo)>>;
    using IntersectionSet = IntersectionEntitySet<GES, AES>;
    const auto gridEntitySet = std::make_shared<GES>(leafGridView);
    const auto averagingEntitySet = std::make_shared<AES>(geo);
    IntersectionSet iset;
    iset.build(gridEntitySet, averagingEntitySet);

    Scalar phiAvg(0.0);
    GlobalPosition vAvg(0.0);
    for(const auto& is : intersections(iset))
    {
        const auto geo = is.geometry();
        const auto volume = geo.volume();
        const auto eIdx = gridEntitySet->index(is.domainEntity());

        phiAvg += volume;
        // Here we apply a piecewise-constant interpolation (because of staggered scheme)
        vAvg += velocity[eIdx]*volume;
    }
    phiAvg /= geo.volume();
    vAvg /= geo.volume();

    Scalar KAvg = vAvg[0]*1e-12;
    std::cout << "Permeability by volume averaging is: " << KAvg << " m^2" << std::endl;
    std::cout << "Porosity by volume averaging is: " << phiAvg << std::endl;

    if (Dune::FloatCmp::ne(KAvg, KRef, 1e-5))
        DUNE_THROW(Dune::Exception, "Permeability (by volume averaging) " << KAvg << " doesn't match reference " << KRef << " by " << KAvg-KRef);
    if (Dune::FloatCmp::ne(phiAvg, phiRef, 1e-5))
        DUNE_THROW(Dune::Exception, "Porosity (by volume averaging)" << phiAvg << " doesn't match reference " << phiRef << " by " << phiAvg-phiRef);

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
