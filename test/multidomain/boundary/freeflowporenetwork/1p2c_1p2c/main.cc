// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoundaryTests
 * \brief A test problem for the coupled Free-Flow/PNM problem (1p).
 */

#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dumux/assembly/diffmethod.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/freeflow/navierstokes/fluxoveraxisalignedsurface.hh>
#include <dumux/freeflow/navierstokes/velocityoutput.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/grid/porenetwork/gridmanager.hh>
#include <dumux/io/vtk/intersectionwriter.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/porenetwork/common/boundaryflux.hh>
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
    using FreeFlowMomentumTypeTag = Properties::TTag::FreeFlowOnePMomentum;
    using FreeFlowMassTypeTag = Properties::TTag::FreeFlowOnePMass;
    using PoreNetworkTypeTag = Properties::TTag::PNMOnePModel;

    using PNMGridManager = Dumux::PoreNetwork::GridManager<2>;
    PNMGridManager pnmGridManager;
    pnmGridManager.init("PNM");

    using FreeFlowGridManager = Dumux::PoreNetwork::SnappyGridManager<2, PNMGridManager>;
    FreeFlowGridManager freeflowGridManager;
    freeflowGridManager.init(pnmGridManager.grid(), *(pnmGridManager.getGridData()), "FreeFlow");
    const auto data = freeflowGridManager.getGridConstructionData();
    const auto auxiliaryPositions = data.interfacePositions[0/*dimIdx*/].value();

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
    auto lowDimspatialParams = std::make_shared<PoreNetworkSpatialParams>(pnmGridGeometry);
    using PNMProblem = GetPropType<PoreNetworkTypeTag, Properties::Problem>;
    auto pnmProblem = std::make_shared<PNMProblem>(pnmGridGeometry, lowDimspatialParams, couplingManager);

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

    // the grid variables
    using FreeFlowMomentumGridVariables = GetPropType<FreeFlowMomentumTypeTag, Properties::GridVariables>;
    auto freeFlowMomentumGridVariables = std::make_shared<FreeFlowMomentumGridVariables>(freeFlowMomentumProblem, freeFlowMomentumGridGeometry);
    using FreeFlowMassGridVariables = GetPropType<FreeFlowMassTypeTag, Properties::GridVariables>;
    auto freeFlowMassGridVariables = std::make_shared<FreeFlowMassGridVariables>(freeFlowMassProblem, freeFlowMassGridGeometry);

    using PoreNetworkgridVariables = GetPropType<PoreNetworkTypeTag, Properties::GridVariables>;
    auto pnmGridVariables = std::make_shared<PoreNetworkgridVariables>(pnmProblem, pnmGridGeometry);

    couplingManager->init(freeFlowMomentumProblem, freeFlowMassProblem, pnmProblem,
                          std::make_tuple(freeFlowMomentumGridVariables, freeFlowMassGridVariables, pnmGridVariables),
                          sol);

    freeFlowMomentumGridVariables->init(sol[freeFlowMomentumIndex]);
    freeFlowMassGridVariables->init(sol[freeFlowMassIndex]);
    pnmGridVariables->init(sol[poreNetworkIndex]);

    // initialize the vtk output module
    VtkOutputModule vtkWriterFF(*freeFlowMassGridVariables, sol[freeFlowMassIndex], freeFlowMassProblem->name());
    GetPropType<FreeFlowMassTypeTag, Properties::IOFields>::initOutputModule(vtkWriterFF); // Add model specific output fields
    vtkWriterFF.addVelocityOutput(std::make_shared<NavierStokesVelocityOutput<FreeFlowMassGridVariables>>());

    using PoreNetworkOutputModule = PoreNetwork::VtkOutputModule<PoreNetworkgridVariables, GetPropType<PoreNetworkTypeTag, Properties::FluxVariables>, std::decay_t<decltype(sol[poreNetworkIndex])>>;
    PoreNetworkOutputModule pmVtkWriter(*pnmGridVariables, sol[poreNetworkIndex], pnmProblem->name());
    GetPropType<PoreNetworkTypeTag, Properties::IOFields>::initOutputModule(pmVtkWriter);

    // the assembler for a stationary problem
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(freeFlowMomentumProblem, freeFlowMassProblem, pnmProblem),
                                                 std::make_tuple(freeFlowMomentumGridGeometry,
                                                                 freeFlowMassGridGeometry,
                                                                 pnmGridGeometry),
                                                 std::make_tuple(freeFlowMomentumGridVariables,
                                                                 freeFlowMassGridVariables,
                                                                 pnmGridVariables),
                                                 couplingManager);

    // the linear solver
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();

    // write vtk output
    vtkWriterFF.write(0.0);
    pmVtkWriter.write(0.0);

    // the non-linear solver
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    // solve the non-linear system
    nonLinearSolver.solve(sol);

    // write vtk output
    vtkWriterFF.write(1.0);
    pmVtkWriter.write(1.0);

    std::vector<Dune::FieldVector<double, 2>> faceVelocity(sol[freeFlowMomentumIndex].size());
    updateVelocities(faceVelocity, *freeFlowMomentumGridGeometry, *freeFlowMomentumGridVariables, sol[freeFlowMomentumIndex]);

    ConformingIntersectionWriter faceVtk(freeFlowMomentumGridGeometry->gridView());
    faceVtk.addField(faceVelocity, "velocityVector");
    faceVtk.write("face_velocities");

    const auto pnmBoundaryFlux = Dumux::PoreNetwork::BoundaryFlux(*pnmGridVariables, assembler->localResidual(poreNetworkIndex), sol[poreNetworkIndex]);
    auto fluxOverSurface = FluxOverAxisAlignedSurface(*freeFlowMassGridVariables, sol[freeFlowMassIndex], assembler->localResidual(freeFlowMassIndex));
    using GlobalPosition = typename FreeFlowMassGridGeometry::SubControlVolume::GlobalPosition;

    if (getParam<bool>("Problem.SingleThroatTest", true))
    {
        const auto lowerLeft = GlobalPosition{freeFlowMassGridGeometry->bBoxMin()[0], freeFlowMassGridGeometry->bBoxMax()[1]};
        const auto upperRight = freeFlowMassGridGeometry->bBoxMax();
        fluxOverSurface.addAxisAlignedSurface("outlet", lowerLeft, upperRight);
        fluxOverSurface.calculateAllFluxes();

        std::cout << "PNM boundary flux at bottom pore is " << pnmBoundaryFlux.getFlux("min", 1, false) << std::endl;
        std::cout << "PNM boundary flux at top pore is " << pnmBoundaryFlux.getFlux("max", 1, false) << std::endl;
        std::cout << "FreeFlow outflow " << fluxOverSurface.flux("outlet") << std::endl;
    }
    else
    {
        int c = 0;
        for (int i = 0; i < auxiliaryPositions.size() - 1; i += 2)
        {
            const auto lowerLeft = GlobalPosition{auxiliaryPositions[i], pnmGridGeometry->bBoxMax()[1]};
            const auto upperRight = GlobalPosition{auxiliaryPositions[i+1], pnmGridGeometry->bBoxMax()[1]};
            fluxOverSurface.addAxisAlignedSurface("couplingPore_" + std::to_string(c++), lowerLeft, upperRight);
        }

        fluxOverSurface.addAxisAlignedPlane("afterFirstPore", 0.25*(pnmGridGeometry->bBoxMin() + pnmGridGeometry->bBoxMax()), 0/*normalAxis*/);
        fluxOverSurface.addAxisAlignedPlane("center", 0.5*(pnmGridGeometry->bBoxMin() + pnmGridGeometry->bBoxMax()), 0/*normalAxis*/);
        fluxOverSurface.addAxisAlignedPlane("outlet", freeFlowMassGridGeometry->bBoxMax(), 0/*normalAxis*/);

        fluxOverSurface.calculateAllFluxes();

        std::cout << "PNM boundary flux at top pores is " << pnmBoundaryFlux.getFlux("max", 1, false) << std::endl;

        fluxOverSurface.printAllFluxes();
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
} // end main
