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
/*!
 * \file
 * \ingroup FacetTests
 * \brief Test for the one-phase facet coupling model.
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
    using BulkFVGridGeometry = Dumux::GetPropType<BulkTypeTag, Dumux::Properties::GridGeometry>;
    using LowDimFVGridGeometry = Dumux::GetPropType<LowDimTypeTag, Dumux::Properties::GridGeometry>;
public:
    using MDTraits = Dumux::MultiDomainTraits<BulkTypeTag, LowDimTypeTag>;
    using CouplingMapper = Dumux::FacetCouplingMapper<BulkFVGridGeometry, LowDimFVGridGeometry>;
    using CouplingManager = Dumux::FacetCouplingManager<MDTraits, CouplingMapper>;
};

// set the coupling manager property in the sub-problems for both box and tpfa
namespace Dumux {
namespace Properties {

// set cm property for the box test
using BoxTraits = TestTraits<Properties::TTag::OnePBulkBox, Properties::TTag::OnePLowDimBox>;
using BoxTracerTraits = TestTraits<Properties::TTag::TracerBulkBox, Properties::TTag::TracerLowDimBox>;
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePBulkBox> { using type = typename BoxTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePLowDimBox> { using type = typename BoxTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::TracerBulkBox> { using type = typename BoxTracerTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::TracerLowDimBox> { using type = typename BoxTracerTraits::CouplingManager; };

// set cm property for the tpfa test
using TpfaTraits = TestTraits<Properties::TTag::OnePBulkTpfa, Properties::TTag::OnePLowDimTpfa>;
using TpfaTracerTraits = TestTraits<Properties::TTag::TracerBulkTpfa, Properties::TTag::TracerLowDimTpfa>;
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePBulkTpfa> { using type = typename TpfaTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePLowDimTpfa> { using type = typename TpfaTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::TracerBulkTpfa> { using type = typename TpfaTracerTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::TracerLowDimTpfa> { using type = typename TpfaTracerTraits::CouplingManager; };

// set cm property for the mpfa test
using MpfaTraits = TestTraits<Properties::TTag::OnePBulkMpfa, Properties::TTag::OnePLowDimMpfa>;
using MpfaTracerTraits = TestTraits<Properties::TTag::TracerBulkMpfa, Properties::TTag::TracerLowDimMpfa>;
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePBulkMpfa> { using type = typename MpfaTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePLowDimMpfa> { using type = typename MpfaTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::TracerBulkMpfa> { using type = typename MpfaTracerTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::TracerLowDimMpfa> { using type = typename MpfaTracerTraits::CouplingManager; };

} // end namespace Properties
} // end namespace Dumux

/*!
 * \brief Updates the finite volume grid geometry for the box scheme.
 *
 * This is necessary as the finite volume grid geometry for the box scheme with
 * facet coupling requires additional data for the update. The reason is that we
 * have to create additional faces on interior boundaries, which wouldn't be
 * created in the standard scheme.
 */
template< class GridGeometry,
          class GridManager,
          class LowDimGridView,
          std::enable_if_t<GridGeometry::discMethod == Dumux::DiscretizationMethod::box, int> = 0 >
void updateBulkFVGridGeometry(GridGeometry& gridGeometry,
                              const GridManager& gridManager,
                              const LowDimGridView& lowDimGridView)
{
    using BulkFacetGridAdapter = Dumux::CodimOneGridAdapter<typename GridManager::Embeddings>;
    BulkFacetGridAdapter facetGridAdapter(gridManager.getEmbeddings());
    gridGeometry.update(lowDimGridView, facetGridAdapter);
}

/*!
 * \brief Updates the finite volume grid geometry for the cell centered schemes.
 */
 template< class GridGeometry,
          class GridManager,
          class LowDimGridView,
          std::enable_if_t<GridGeometry::discMethod != Dumux::DiscretizationMethod::box, int> = 0 >
void updateBulkFVGridGeometry(GridGeometry& gridGeometry,
                              const GridManager& gridManager,
                              const LowDimGridView& lowDimGridView)
{
    gridGeometry.update();
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
    static constexpr bool isBox = GV::GridGeometry::discMethod == Dumux::DiscretizationMethod::box;

    // resize depending on the scheme
    if (!isBox) volumeFluxes.assign(gridGeometry.numScvf(), {0.0});
    else volumeFluxes.assign(gridGeometry.gridView().size(0), {0.0});

    auto upwindTerm = [](const auto& volVars) { return volVars.mobility(0); };
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        const auto eIdx = gridGeometry.elementMapper().index(element);

        // bind local views
        couplingManager.bindCouplingContext(domainId, element, assembler);
        auto fvGeometry = localView(gridGeometry);
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

// main program
int main(int argc, char** argv)
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

        // intialize the vtk output module
        const auto bulkDM = BulkFVGridGeometry::discMethod == DiscretizationMethod::box ? Dune::VTK::nonconforming : Dune::VTK::conforming;
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
    using TracerTestTraits = TestTraits<BulkTracerTypeTag, LowDimTracerTypeTag>;
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

    // intialize the vtk output modules
    const auto bulkDM = BulkFVGridGeometry::discMethod == DiscretizationMethod::box ? Dune::VTK::nonconforming : Dune::VTK::conforming;
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
