// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
/*!
 * \file
 *
 * \brief test for the one-phase facet coupling model
 */
#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/geometry/quadraturerules.hh>

#include "problem_1p_bulk.hh"
#include "problem_1p_lowdim.hh"

#include "problem_tracer_bulk.hh"
#include "problem_tracer_lowdim.hh"

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>

#include <dumux/assembly/diffmethod.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/evalgradients.hh>
#include <dumux/linear/seqsolverbackend.hh>

#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>

#include <dumux/multidomain/facet/gridmanager.hh>
#include <dumux/multidomain/facet/couplingmapper.hh>
#include <dumux/multidomain/facet/couplingmanager.hh>
#include <dumux/multidomain/facet/codimonegridadapter.hh>

#include <dumux/io/vtkoutputmodule.hh>

// obtain/define some types to be used below in the property definitions and in main
template< class BulkTypeTag, class LowDimTypeTag >
class TestTraits
{
    using BulkFVGridGeometry = typename GET_PROP_TYPE(BulkTypeTag, FVGridGeometry);
    using LowDimFVGridGeometry = typename GET_PROP_TYPE(LowDimTypeTag, FVGridGeometry);
public:
    using MDTraits = Dumux::MultiDomainTraits<BulkTypeTag, LowDimTypeTag>;
    using CouplingMapper = Dumux::FacetCouplingMapper<BulkFVGridGeometry, LowDimFVGridGeometry>;
    using CouplingManager = Dumux::FacetCouplingManager<MDTraits, CouplingMapper>;
};

// set the coupling manager property in the sub-problems for both box and tpfa
namespace Dumux {
namespace Properties {

// set cm property for the box test
using BoxTraits = TestTraits<TTAG(OnePBulkBox), TTAG(OnePLowDimBox)>;
using BoxTracerTraits = TestTraits<TTAG(TracerBulkBox), TTAG(TracerLowDimBox)>;
SET_TYPE_PROP(OnePBulkBox, CouplingManager, typename BoxTraits::CouplingManager);
SET_TYPE_PROP(OnePLowDimBox, CouplingManager, typename BoxTraits::CouplingManager);
SET_TYPE_PROP(TracerBulkBox, CouplingManager, typename BoxTracerTraits::CouplingManager);
SET_TYPE_PROP(TracerLowDimBox, CouplingManager, typename BoxTracerTraits::CouplingManager);

// set cm property for the tpfa test
using TpfaTraits = TestTraits<TTAG(OnePBulkTpfa), TTAG(OnePLowDimTpfa)>;
using TpfaTracerTraits = TestTraits<TTAG(TracerBulkTpfa), TTAG(TracerLowDimTpfa)>;
SET_TYPE_PROP(OnePBulkTpfa, CouplingManager, typename TpfaTraits::CouplingManager);
SET_TYPE_PROP(OnePLowDimTpfa, CouplingManager, typename TpfaTraits::CouplingManager);
SET_TYPE_PROP(TracerBulkTpfa, CouplingManager, typename TpfaTracerTraits::CouplingManager);
SET_TYPE_PROP(TracerLowDimTpfa, CouplingManager, typename TpfaTracerTraits::CouplingManager);

} // end namespace Properties
} // end namespace Dumux

// updates the finite volume grid geometry. This is necessary as the finite volume
// grid geometry for the box scheme with facet coupling requires additional data for
// the update. The reason is that we have to create additional faces on interior
// boundaries, which wouldn't be created in the standard scheme.
template< class FVGridGeometry,
          class GridManager,
          class LowDimGridView,
          std::enable_if_t<FVGridGeometry::discMethod == Dumux::DiscretizationMethod::box, int> = 0 >
void updateBulkFVGridGeometry(FVGridGeometry& fvGridGeometry,
                              const GridManager& gridManager,
                              const LowDimGridView& lowDimGridView)
{
    using BulkFacetGridAdapter = Dumux::CodimOneGridAdapter<typename GridManager::Embeddings>;
    BulkFacetGridAdapter facetGridAdapter(gridManager.getEmbeddings());
    fvGridGeometry.update(lowDimGridView, facetGridAdapter);
}

// specialization for cell-centered schemes
template< class FVGridGeometry,
          class GridManager,
          class LowDimGridView,
          std::enable_if_t<FVGridGeometry::discMethod != Dumux::DiscretizationMethod::box, int> = 0 >
void updateBulkFVGridGeometry(FVGridGeometry& fvGridGeometry,
                              const GridManager& gridManager,
                              const LowDimGridView& lowDimGridView)
{
    fvGridGeometry.update();
}

//! computes the volume fluxes on all scvfs for a sub-domain
template<class FV, class Storage, class CM, class Assembler, class Prob, class GV, class Sol, std::size_t id>
void computeVolumeFluxes(Storage& volumeFluxes,
                         CM& couplingManager,
                         const Assembler& assembler,
                         const Prob& problem,
                         const typename GV::GridGeometry& fvGridGeometry,
                         const GV& gridVariables,
                         const Sol& sol,
                         Dune::index_constant<id> domainId)
{
    static constexpr bool isBox = GV::GridGeometry::discMethod == Dumux::DiscretizationMethod::box;

    // resize depending on the scheme
    if (!isBox) volumeFluxes.assign(fvGridGeometry.numScvf(), {0.0});
    else volumeFluxes.assign(fvGridGeometry.gridView().size(0), {0.0});

    auto upwindTerm = [](const auto& volVars) { return volVars.mobility(0); };
    for (const auto& element : elements(fvGridGeometry.gridView()))
    {
        const auto eIdx = fvGridGeometry.elementMapper().index(element);

        // bind local views
        couplingManager.bindCouplingContext(domainId, element, assembler);
        auto fvGeometry = localView(fvGridGeometry);
        auto elemVolVars = localView(gridVariables.curGridVolVars());
        auto elemFluxVars = localView(gridVariables.gridFluxVarsCache());
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
                const auto elemSol = elementSolution(element, sol, fvGridGeometry);
                const auto gradP = evalGradients(element,
                                                 element.geometry(),
                                                 fvGridGeometry,
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

// main program
int main(int argc, char** argv) try
{
    using namespace Dumux;

    //////////////////////////////////////////////////////
    //////////////////////////////////////////////////////

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // initialize parameter tree
    Parameters::init(argc, argv);

    //////////////////////////////////////////////////////
    // try to create the grids (from the given grid file)
    //////////////////////////////////////////////////////
    using BulkOnePTypeTag = TTAG(ONEPBULKTYPETAG);
    using LowDimOnePTypeTag = TTAG(ONEPLOWDIMTYPETAG);
    using BulkGrid = typename GET_PROP_TYPE(BulkOnePTypeTag, Grid);
    using LowDimGrid = typename GET_PROP_TYPE(LowDimOnePTypeTag, Grid);

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
    using BulkFVGridGeometry = typename GET_PROP_TYPE(BulkOnePTypeTag, FVGridGeometry);
    using LowDimFVGridGeometry = typename GET_PROP_TYPE(LowDimOnePTypeTag, FVGridGeometry);
    auto bulkFvGridGeometry = std::make_shared<BulkFVGridGeometry>(bulkGridView);
    auto lowDimFvGridGeometry = std::make_shared<LowDimFVGridGeometry>(lowDimGridView);
    updateBulkFVGridGeometry(*bulkFvGridGeometry, gridManager, lowDimGridView);
    lowDimFvGridGeometry->update();

    // the coupling mapper
    using OnePTestTraits = TestTraits<BulkOnePTypeTag, LowDimOnePTypeTag>;
    auto couplingMapper = std::make_shared<typename OnePTestTraits::CouplingMapper>();
    couplingMapper->update(*bulkFvGridGeometry, *lowDimFvGridGeometry, gridManager.getEmbeddings());

    // get the ids of bulk/facet domain & linear system types
    // they are the same types for onep/tracer - obtain only once here
    using OnePMDTraits = typename OnePTestTraits::MDTraits;
    using SolutionVector = typename OnePMDTraits::SolutionVector;
    using JacobianMatrix = typename OnePMDTraits::JacobianMatrix;
    static const auto bulkId = typename OnePMDTraits::template DomainIdx<0>();
    static const auto lowDimId = typename OnePMDTraits::template DomainIdx<1>();

    ////////////////////////////////////////////////////////////
    // run stationary, simgle-phase problem on this grid
    ////////////////////////////////////////////////////////////

    // put the following code in brackets such that memory is released afterwards
    {
        // instantiate coupling manager
        using CouplingManager = typename OnePTestTraits::CouplingManager;
        auto couplingManager = std::make_shared<CouplingManager>();

        // the problems (boundary conditions)
        using BulkProblem = typename GET_PROP_TYPE(BulkOnePTypeTag, Problem);
        using LowDimProblem = typename GET_PROP_TYPE(LowDimOnePTypeTag, Problem);
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
        using BulkGridVariables = typename GET_PROP_TYPE(BulkOnePTypeTag, GridVariables);
        using LowDimGridVariables = typename GET_PROP_TYPE(LowDimOnePTypeTag, GridVariables);
        auto bulkGridVariables = std::make_shared<BulkGridVariables>(bulkProblem, bulkFvGridGeometry);
        auto lowDimGridVariables = std::make_shared<LowDimGridVariables>(lowDimProblem, lowDimFvGridGeometry);
        bulkGridVariables->init(x[bulkId]);
        lowDimGridVariables->init(x[lowDimId]);

        // intialize the vtk output module
        const auto bulkDM = BulkFVGridGeometry::discMethod == DiscretizationMethod::box ? Dune::VTK::nonconforming : Dune::VTK::conforming;
        using BulkSolutionVector = std::decay_t<decltype(x[bulkId])>;
        using LowDimSolutionVector = std::decay_t<decltype(x[lowDimId])>;
        VtkOutputModule<BulkGridVariables, BulkSolutionVector> bulkVtkWriter(*bulkGridVariables, x[bulkId], bulkProblem->name(), "Bulk.OneP", bulkDM);
        VtkOutputModule<LowDimGridVariables, LowDimSolutionVector> lowDimVtkWriter(*lowDimGridVariables, x[lowDimId], lowDimProblem->name(), "LowDim.OneP");

        // Add model specific output fields
        using BulkIOFields = typename GET_PROP_TYPE(BulkOnePTypeTag, IOFields);
        using LowDimIOFields = typename GET_PROP_TYPE(LowDimOnePTypeTag, IOFields);
        BulkIOFields::initOutputModule(bulkVtkWriter);
        LowDimIOFields::initOutputModule(lowDimVtkWriter);

        // the assembler
        using Assembler = MultiDomainFVAssembler<OnePMDTraits, CouplingManager, DiffMethod::numeric, /*implicit?*/true>;
        auto assembler = std::make_shared<Assembler>( std::make_tuple(bulkProblem, lowDimProblem),
                                                      std::make_tuple(bulkFvGridGeometry, lowDimFvGridGeometry),
                                                      std::make_tuple(bulkGridVariables, lowDimGridVariables),
                                                      couplingManager);

        // the linear solver
        using LinearSolver = UMFPackBackend;
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
        using BulkFluxVariables = typename GET_PROP_TYPE(BulkOnePTypeTag, FluxVariables);
        using LowDimFluxVariables = typename GET_PROP_TYPE(LowDimOnePTypeTag, FluxVariables);
        computeVolumeFluxes<BulkFluxVariables>(bulkVolumeFluxes, *couplingManager, *assembler,
                                               *bulkProblem, *bulkFvGridGeometry, *bulkGridVariables, x[bulkId], bulkId);
        computeVolumeFluxes<LowDimFluxVariables>(lowDimVolumeFluxes, *couplingManager, *assembler,
                                                 *lowDimProblem, *lowDimFvGridGeometry, *lowDimGridVariables, x[lowDimId], lowDimId);
    }

    ////////////////////////////////////////////////////////////////////////////
    // run instationary tracer problem on this grid with the precomputed fluxes
    ////////////////////////////////////////////////////////////////////////////

    //! the problem (initial and boundary conditions)
    using BulkTracerTypeTag = TTAG(TRACERBULKTYPETAG);
    using LowDimTracerTypeTag = TTAG(TRACERLOWDIMTYPETAG);

    // instantiate coupling manager
    using TracerTestTraits = TestTraits<BulkTracerTypeTag, LowDimTracerTypeTag>;
    using CouplingManager = typename TracerTestTraits::CouplingManager;
    auto couplingManager = std::make_shared<CouplingManager>();

    // instantiate the tracer problems reusing the fv grid geometries
    using TracerBulkProblem = typename GET_PROP_TYPE(BulkTracerTypeTag, Problem);
    using TracerLowDimProblem = typename GET_PROP_TYPE(LowDimTracerTypeTag, Problem);
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
    using BulkGridVariables = typename GET_PROP_TYPE(BulkTracerTypeTag, GridVariables);
    using LowDimGridVariables = typename GET_PROP_TYPE(LowDimTracerTypeTag, GridVariables);
    auto bulkGridVariables = std::make_shared<BulkGridVariables>(bulkProblem, bulkFvGridGeometry);
    auto lowDimGridVariables = std::make_shared<LowDimGridVariables>(lowDimProblem, lowDimFvGridGeometry);
    bulkGridVariables->init(x[bulkId]);
    lowDimGridVariables->init(x[lowDimId]);

    // intialize the vtk output modules
    const auto bulkDM = BulkFVGridGeometry::discMethod == DiscretizationMethod::box ? Dune::VTK::nonconforming : Dune::VTK::conforming;
    using BulkSolutionVector = std::decay_t<decltype(x[bulkId])>;
    using LowDimSolutionVector = std::decay_t<decltype(x[lowDimId])>;
    VtkOutputModule<BulkGridVariables, BulkSolutionVector> bulkVtkWriter(*bulkGridVariables, x[bulkId], bulkProblem->name(), "Bulk.Tracer", bulkDM);
    VtkOutputModule<LowDimGridVariables, LowDimSolutionVector> lowDimVtkWriter(*lowDimGridVariables, x[lowDimId], lowDimProblem->name(), "LowDim.Tracer");

    // Add model specific output fields
    using BulkIOFields = typename GET_PROP_TYPE(BulkTracerTypeTag, IOFields);
    using LowDimIOFields = typename GET_PROP_TYPE(LowDimTracerTypeTag, IOFields);
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
                                                  couplingManager, timeLoop);

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
        assembler->setPreviousSolution(xOld);
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
catch (Dumux::ParameterException &e)
{
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
}
catch (Dune::DGFException & e)
{
    std::cerr << "DGF exception thrown (" << e <<
                 "). Most likely, the DGF file name is wrong "
                 "or the DGF file is corrupted, "
                 "e.g. missing hash at end of file or wrong number (dimensions) of entries."
                 << " ---> Abort!" << std::endl;
    return 2;
}
catch (Dune::Exception &e)
{
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
}
catch (...)
{
    std::cerr << "Unknown exception thrown! ---> Abort!" << std::endl;
    return 4;
}
