// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FacetTests
 * \brief Test for the one-phase facet coupling model.
 */

#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/evalgradients.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>

#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>

#include <dumux/multidomain/facet/gridmanager.hh>
#include <dumux/multidomain/facet/codimonegridadapter.hh>

#include <dumux/io/vtkoutputmodule.hh>

#include "properties.hh"

namespace Dumux {

/*!
 * \brief Constructs the finite volume grid geometry.
 */
template< class BulkGridGeometry,
          class GridManager,
          class BulkGridView,
          class LowDimGridView >
auto makeBulkFVGridGeometry(const GridManager& gridManager,
                            const BulkGridView& bulkGridView,
                            const LowDimGridView& lowDimGridView)
{
    /*!
    * The finite volume grid geometry for the box scheme with facet coupling
    * requires additional data for the constructor. The reason is that
    * we have to create additional faces on interior boundaries, which are not
    * created in the standard scheme.
    */
    if constexpr (BulkGridGeometry::discMethod == Dumux::DiscretizationMethods::box)
    {
        using BulkFacetGridAdapter = Dumux::CodimOneGridAdapter<typename GridManager::Embeddings>;
        BulkFacetGridAdapter facetGridAdapter(gridManager.getEmbeddings());
        return std::make_shared<BulkGridGeometry>(bulkGridView, lowDimGridView, facetGridAdapter);
    }
    else
    {
        return std::make_shared<BulkGridGeometry>(bulkGridView);
    }
}

//! Computes the volume fluxes on all scvfs for a sub-domain.
template<class FV, class Storage, class CM, class Assembler, class Prob, class GV, class Sol, std::size_t id>
void computeVolumeFluxes(Storage& volumeFluxes,
                         CM& couplingManager,
                         const Assembler& assembler,
                         const Prob& problem,
                         const typename GV::GridGeometry& gridGeometry,
                         const GV& gridVariables,
                         const Sol& sol,
                         Dune::index_constant<id> domainId)
{
    static constexpr bool isBox = GV::GridGeometry::discMethod == Dumux::DiscretizationMethods::box;

    // resize depending on the scheme
    if (!isBox) volumeFluxes.assign(gridGeometry.numScvf(), {0.0});
    else volumeFluxes.assign(gridGeometry.gridView().size(0), {0.0});

    auto upwindTerm = [](const auto& volVars) { return volVars.mobility(0); };
    auto fvGeometry = localView(gridGeometry);
    auto elemVolVars = localView(gridVariables.curGridVolVars());
    auto elemFluxVars = localView(gridVariables.gridFluxVarsCache());
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        const auto eIdx = gridGeometry.elementMapper().index(element);

        // bind local views
        couplingManager.bindCouplingContext(domainId, element, assembler);
        fvGeometry.bind(element);
        elemVolVars.bind(element, fvGeometry, sol);
        elemFluxVars.bind(element, fvGeometry, elemVolVars);

        if (isBox)
            volumeFluxes[eIdx].resize(fvGeometry.numScvf(), 0.0);

        for (const auto& scvf : scvfs(fvGeometry))
        {
            typename GV::Scalar flux = 0.0;
            FV fluxVars;
            fluxVars.init(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVars);

            if (!scvf.boundary())
                flux = fluxVars.advectiveFlux(0, upwindTerm);
            else if (!isBox && problem.boundaryTypes(element, scvf).hasOnlyDirichlet())
                flux = fluxVars.advectiveFlux(0, upwindTerm);
            else if (isBox
                     && problem.boundaryTypes(element, fvGeometry.scv(scvf.insideScvIdx())).hasOnlyDirichlet())
            {
                // reconstruct flux
                const auto elemSol = elementSolution(element, sol, gridGeometry);
                const auto gradP = evalGradients(element,
                                                 element.geometry(),
                                                 gridGeometry,
                                                 elemSol,
                                                 scvf.ipGlobal())[0];
                const auto& insideVolVars = elemVolVars[fvGeometry.scv(scvf.insideScvIdx())];
                const auto k = insideVolVars.permeability();
                flux = - insideVolVars.mobility(0)*k*(gradP*scvf.unitOuterNormal())*scvf.area();
            }

            if (isBox) volumeFluxes[eIdx][scvf.index()] = flux;
            else volumeFluxes[scvf.index()][0] = flux;
        }
    }
}

} // end namespace Dumux

// main program
int main(int argc, char** argv)
{
    using namespace Dumux;

    //////////////////////////////////////////////////////
    //////////////////////////////////////////////////////

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // initialize parameter tree
    Parameters::init(argc, argv);

    //////////////////////////////////////////////////////
    // try to create the grids (from the given grid file)
    //////////////////////////////////////////////////////
    using BulkOnePTypeTag = Properties::TTag::ONEPBULKTYPETAG;
    using LowDimOnePTypeTag = Properties::TTag::ONEPLOWDIMTYPETAG;
    using BulkGrid = GetPropType<BulkOnePTypeTag, Properties::Grid>;
    using LowDimGrid = GetPropType<LowDimOnePTypeTag, Properties::Grid>;

    using GridManager = FacetCouplingGridManager<BulkGrid, LowDimGrid>;
    GridManager gridManager;
    gridManager.init();
    gridManager.loadBalance();

    // we compute on the leaf grid views
    const auto& bulkGridView = gridManager.template grid<0>().leafGridView();
    const auto& lowDimGridView = gridManager.template grid<1>().leafGridView();

    // containers to store the volume fluxes (will be passed to tracer problem)
    std::vector< std::vector<double> > bulkVolumeFluxes;
    std::vector< std::vector<double> > lowDimVolumeFluxes;

    // create the finite volume grid geometries
    using BulkFVGridGeometry = GetPropType<BulkOnePTypeTag, Properties::GridGeometry>;
    using LowDimFVGridGeometry = GetPropType<LowDimOnePTypeTag, Properties::GridGeometry>;
    auto lowDimFvGridGeometry = std::make_shared<LowDimFVGridGeometry>(lowDimGridView);
    auto bulkFvGridGeometry  = makeBulkFVGridGeometry<BulkFVGridGeometry>(gridManager, bulkGridView, lowDimGridView);

    // the coupling mapper
    using OnePTestTraits = Properties::TestTraits<BulkOnePTypeTag, LowDimOnePTypeTag>;
    auto couplingMapper = std::make_shared<typename OnePTestTraits::CouplingMapper>();
    couplingMapper->update(*bulkFvGridGeometry, *lowDimFvGridGeometry, gridManager.getEmbeddings());

    // get the ids of bulk/facet domain & linear system types
    // they are the same types for onep/tracer - obtain only once here
    using OnePMDTraits = typename OnePTestTraits::MDTraits;
    using SolutionVector = typename OnePMDTraits::SolutionVector;
    using JacobianMatrix = typename OnePMDTraits::JacobianMatrix;
    static const auto bulkId = typename OnePMDTraits::template SubDomain<0>::Index();
    static const auto lowDimId = typename OnePMDTraits::template SubDomain<1>::Index();

    ////////////////////////////////////////////////////////////
    // run stationary, simgle-phase problem on this grid
    ////////////////////////////////////////////////////////////

    // put the following code in brackets such that memory is released afterwards
    {
        // instantiate coupling manager
        using CouplingManager = typename OnePTestTraits::CouplingManager;
        auto couplingManager = std::make_shared<CouplingManager>();

        // the problems (boundary conditions)
        using BulkProblem = GetPropType<BulkOnePTypeTag, Properties::Problem>;
        using LowDimProblem = GetPropType<LowDimOnePTypeTag, Properties::Problem>;
        auto bulkSpatialParams = std::make_shared<typename BulkProblem::SpatialParams>(bulkFvGridGeometry, "Bulk.OneP");
        auto bulkProblem = std::make_shared<BulkProblem>(bulkFvGridGeometry, bulkSpatialParams, couplingManager, "Bulk.OneP");
        auto lowDimSpatialParams = std::make_shared<typename LowDimProblem::SpatialParams>(lowDimFvGridGeometry, "LowDim.OneP");
        auto lowDimProblem = std::make_shared<LowDimProblem>(lowDimFvGridGeometry, lowDimSpatialParams, couplingManager, "LowDim.OneP");

        // the solution vector
        SolutionVector x;
        x[bulkId].resize(bulkFvGridGeometry->numDofs());
        x[lowDimId].resize(lowDimFvGridGeometry->numDofs());
        bulkProblem->applyInitialSolution(x[bulkId]);
        lowDimProblem->applyInitialSolution(x[lowDimId]);

        // initialize the coupling manager
        couplingManager->init(bulkProblem, lowDimProblem, couplingMapper, x);

        // the grid variables
        using BulkGridVariables = GetPropType<BulkOnePTypeTag, Properties::GridVariables>;
        using LowDimGridVariables = GetPropType<LowDimOnePTypeTag, Properties::GridVariables>;
        auto bulkGridVariables = std::make_shared<BulkGridVariables>(bulkProblem, bulkFvGridGeometry);
        auto lowDimGridVariables = std::make_shared<LowDimGridVariables>(lowDimProblem, lowDimFvGridGeometry);
        bulkGridVariables->init(x[bulkId]);
        lowDimGridVariables->init(x[lowDimId]);

        // initialize the vtk output module
        const auto bulkDM = BulkFVGridGeometry::discMethod == DiscretizationMethods::box ? Dune::VTK::nonconforming : Dune::VTK::conforming;
        using BulkSolutionVector = std::decay_t<decltype(x[bulkId])>;
        using LowDimSolutionVector = std::decay_t<decltype(x[lowDimId])>;
        VtkOutputModule<BulkGridVariables, BulkSolutionVector> bulkVtkWriter(*bulkGridVariables, x[bulkId], bulkProblem->name(), "Bulk.OneP", bulkDM);
        VtkOutputModule<LowDimGridVariables, LowDimSolutionVector> lowDimVtkWriter(*lowDimGridVariables, x[lowDimId], lowDimProblem->name(), "LowDim.OneP");

        // Add model specific output fields
        using BulkIOFields = GetPropType<BulkOnePTypeTag, Properties::IOFields>;
        using LowDimIOFields = GetPropType<LowDimOnePTypeTag, Properties::IOFields>;
        BulkIOFields::initOutputModule(bulkVtkWriter);
        LowDimIOFields::initOutputModule(lowDimVtkWriter);

        // the assembler
        using Assembler = MultiDomainFVAssembler<OnePMDTraits, CouplingManager, DiffMethod::numeric, /*implicit?*/true>;
        auto assembler = std::make_shared<Assembler>( std::make_tuple(bulkProblem, lowDimProblem),
                                                      std::make_tuple(bulkFvGridGeometry, lowDimFvGridGeometry),
                                                      std::make_tuple(bulkGridVariables, lowDimGridVariables),
                                                      couplingManager);

        // the linear solver
        using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
        auto linearSolver = std::make_shared<LinearSolver>();

        // the non-linear solver
        using NewtonSolver = Dumux::MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
        auto newtonSolver = std::make_shared<NewtonSolver>(assembler, linearSolver, couplingManager);

        // linearize & solve
        newtonSolver->solve(x);

        // update grid variables for output
        bulkGridVariables->update(x[bulkId]);
        lowDimGridVariables->update(x[lowDimId]);

        // write vtk output
        bulkVtkWriter.write(1.0);
        lowDimVtkWriter.write(1.0);

        // compute the volume fluxes and store them in the arrays
        using BulkFluxVariables = GetPropType<BulkOnePTypeTag, Properties::FluxVariables>;
        using LowDimFluxVariables = GetPropType<LowDimOnePTypeTag, Properties::FluxVariables>;
        computeVolumeFluxes<BulkFluxVariables>(bulkVolumeFluxes, *couplingManager, *assembler,
                                               *bulkProblem, *bulkFvGridGeometry, *bulkGridVariables, x[bulkId], bulkId);
        computeVolumeFluxes<LowDimFluxVariables>(lowDimVolumeFluxes, *couplingManager, *assembler,
                                                 *lowDimProblem, *lowDimFvGridGeometry, *lowDimGridVariables, x[lowDimId], lowDimId);
    }

    ////////////////////////////////////////////////////////////////////////////
    // run instationary tracer problem on this grid with the precomputed fluxes
    ////////////////////////////////////////////////////////////////////////////

    //! the problem (initial and boundary conditions)
    using BulkTracerTypeTag = Properties::TTag::TRACERBULKTYPETAG;
    using LowDimTracerTypeTag = Properties::TTag::TRACERLOWDIMTYPETAG;

    // instantiate coupling manager
    using TracerTestTraits = Properties::TestTraits<BulkTracerTypeTag, LowDimTracerTypeTag>;
    using CouplingManager = typename TracerTestTraits::CouplingManager;
    auto couplingManager = std::make_shared<CouplingManager>();

    // instantiate the tracer problems reusing the fv grid geometries
    using TracerBulkProblem = GetPropType<BulkTracerTypeTag, Properties::Problem>;
    using TracerLowDimProblem = GetPropType<LowDimTracerTypeTag, Properties::Problem>;
    using TracerBulkSpatialParams = typename TracerBulkProblem::SpatialParams;
    using TracerLowDimSpatialParams = typename TracerLowDimProblem::SpatialParams;
    auto bulkSpatialParams = std::make_shared<TracerBulkSpatialParams>(bulkFvGridGeometry, bulkVolumeFluxes, "Bulk.Tracer");
    auto bulkProblem = std::make_shared<TracerBulkProblem>(bulkFvGridGeometry, bulkSpatialParams, couplingManager, "Bulk.Tracer");
    auto lowDimSpatialParams = std::make_shared<TracerLowDimSpatialParams>(lowDimFvGridGeometry, lowDimVolumeFluxes, "LowDim.Tracer");
    auto lowDimProblem = std::make_shared<TracerLowDimProblem>(lowDimFvGridGeometry, lowDimSpatialParams, couplingManager, "LowDim.Tracer");

    // the solution vectors and system matrix
    SolutionVector x, xOld;
    x[bulkId].resize(bulkFvGridGeometry->numDofs());
    x[lowDimId].resize(lowDimFvGridGeometry->numDofs());
    bulkProblem->applyInitialSolution(x[bulkId]);
    lowDimProblem->applyInitialSolution(x[lowDimId]);
    xOld = x;

    //! the linear system
    auto A = std::make_shared<JacobianMatrix>();
    auto r = std::make_shared<SolutionVector>();

    // initialize the coupling manager
    couplingManager->init(bulkProblem, lowDimProblem, couplingMapper, x);

    // the grid variables
    using BulkGridVariables = GetPropType<BulkTracerTypeTag, Properties::GridVariables>;
    using LowDimGridVariables = GetPropType<LowDimTracerTypeTag, Properties::GridVariables>;
    auto bulkGridVariables = std::make_shared<BulkGridVariables>(bulkProblem, bulkFvGridGeometry);
    auto lowDimGridVariables = std::make_shared<LowDimGridVariables>(lowDimProblem, lowDimFvGridGeometry);
    bulkGridVariables->init(x[bulkId]);
    lowDimGridVariables->init(x[lowDimId]);

    // initialize the vtk output modules
    const auto bulkDM = BulkFVGridGeometry::discMethod == DiscretizationMethods::box ? Dune::VTK::nonconforming : Dune::VTK::conforming;
    using BulkSolutionVector = std::decay_t<decltype(x[bulkId])>;
    using LowDimSolutionVector = std::decay_t<decltype(x[lowDimId])>;
    VtkOutputModule<BulkGridVariables, BulkSolutionVector> bulkVtkWriter(*bulkGridVariables, x[bulkId], bulkProblem->name(), "Bulk.Tracer", bulkDM);
    VtkOutputModule<LowDimGridVariables, LowDimSolutionVector> lowDimVtkWriter(*lowDimGridVariables, x[lowDimId], lowDimProblem->name(), "LowDim.Tracer");

    // Add model specific output fields
    using BulkIOFields = GetPropType<BulkTracerTypeTag, Properties::IOFields>;
    using LowDimIOFields = GetPropType<LowDimTracerTypeTag, Properties::IOFields>;
    BulkIOFields::initOutputModule(bulkVtkWriter);
    LowDimIOFields::initOutputModule(lowDimVtkWriter);

    // write initial solution
    bulkVtkWriter.write(0.0);
    lowDimVtkWriter.write(0.0);

    //! get some time loop parameters
    const auto tEnd = getParam<double>("TimeLoop.TEnd");
    auto dt = getParam<double>("TimeLoop.Dt");

    //! instantiate time loop
    auto timeLoop = std::make_shared<CheckPointTimeLoop<double>>(0.0, dt, tEnd);

    // the assembler
    using Assembler = MultiDomainFVAssembler<typename TracerTestTraits::MDTraits, CouplingManager, DiffMethod::numeric, /*implicit?*/false>;
    auto assembler = std::make_shared<Assembler>( std::make_tuple(bulkProblem, lowDimProblem),
                                                  std::make_tuple(bulkFvGridGeometry, lowDimFvGridGeometry),
                                                  std::make_tuple(bulkGridVariables, lowDimGridVariables),
                                                  couplingManager, timeLoop, xOld);

    // tell the assembler which matrix/residual to use
    assembler->setLinearSystem(A, r);

    // the linear solver
    using LinearSolver = BlockDiagAMGBiCGSTABSolver;
    auto linearSolver = std::make_shared<LinearSolver>();

    //! set some check points for the vtk output
    using std::max;
    timeLoop->setPeriodicCheckPoint( max(tEnd/10.0, dt) );

    //! start the time loop
    timeLoop->start(); do
    {
        // set previous solution for storage evaluations
        Dune::Timer assembleTimer;
        assembler->assembleJacobianAndResidual(x);
        assembleTimer.stop();

        // solve the linear system A(xOld-xNew) = r
        Dune::Timer solveTimer;
        SolutionVector xDelta(x);
        linearSolver->solve(*A, xDelta, *r);
        solveTimer.stop();

        // update solution and grid variables
        Dune::Timer updateTimer;
        x -= xDelta;
        bulkGridVariables->update(x[bulkId]);
        lowDimGridVariables->update(x[lowDimId]);
        updateTimer.stop();

        // statistics
        const auto elapsedTot = assembleTimer.elapsed() + solveTimer.elapsed() + updateTimer.elapsed();
        std::cout << "Assemble/solve/update time: "
                  <<  assembleTimer.elapsed() << "(" << 100*assembleTimer.elapsed()/elapsedTot << "%)/"
                  <<  solveTimer.elapsed() << "(" << 100*solveTimer.elapsed()/elapsedTot << "%)/"
                  <<  updateTimer.elapsed() << "(" << 100*updateTimer.elapsed()/elapsedTot << "%)"
                  <<  std::endl;

        // make the new solution the old solution
        xOld = x;
        bulkGridVariables->advanceTimeStep();
        lowDimGridVariables->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // write vtk output on check points
        if (timeLoop->isCheckPoint())
        {
            bulkVtkWriter.write(timeLoop->time());
            lowDimVtkWriter.write(timeLoop->time());
        }

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt
        timeLoop->setTimeStepSize(dt);

    } while (!timeLoop->finished());

    timeLoop->finalize();

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/false);

    return 0;
}
