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
 * \brief Test for the compositional one-phase facet coupling model.
 */
#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/assembly/diffmethod.hh>
#include <dumux/linear/seqsolverbackend.hh>

#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>

#include <dumux/multidomain/facet/gridmanager.hh>
#include <dumux/multidomain/facet/couplingmapper.hh>
#include <dumux/multidomain/facet/couplingmanager.hh>
#include <dumux/multidomain/facet/codimonegridadapter.hh>

#include <dumux/io/vtkoutputmodule.hh>

#include "properties_bulk.hh"
#include "properties_facet.hh"

using BulkTypeTag = Dumux::Properties::TTag::BULKTYPETAG;
using FacetTypeTag = Dumux::Properties::TTag::FACETTYPETAG;

namespace Dumux {

// obtain/define some types to be used below in the property definitions and in main
class TestTraits
{
    using BulkGridGeometry = GetPropType<BulkTypeTag, Properties::GridGeometry>;
    using FacetGridGeometry = GetPropType<FacetTypeTag, Properties::GridGeometry>;
public:
    using MDTraits = Dumux::MultiDomainTraits<BulkTypeTag, FacetTypeTag>;
    using CouplingMapper = Dumux::FacetCouplingMapper<BulkGridGeometry, FacetGridGeometry>;
    using CouplingManager = Dumux::FacetCouplingManager<MDTraits, CouplingMapper>;
};

// specify coupling manager property in sub-problems
namespace Properties {
template<class TypeTag> struct CouplingManager<TypeTag, BulkTypeTag> { using type = typename TestTraits::CouplingManager; };
template<class TypeTag> struct CouplingManager<TypeTag, FacetTypeTag> { using type = typename TestTraits::CouplingManager; };
} // end namespace Properties
} // end namespace Dumux

/*!
 * \brief Updates the finite volume grid geometry for the box scheme.
 *
 * This is necessary as the finite volume grid geometry for the box scheme with
 * facet coupling requires additional data for the update. The reason is that
 * we have to create additional faces on interior boundaries, which wouldn't be
 * created in the standard scheme.
 */
template< class GridGeometry,
          class GridManager,
          class FacetGridView,
          std::enable_if_t<GridGeometry::discMethod == Dumux::DiscretizationMethod::box, int> = 0 >
void updateBulkFVGridGeometry(GridGeometry& gridGeometry,
                              const GridManager& gridManager,
                              const FacetGridView& facetGridView)
{
    using BulkFacetGridAdapter = Dumux::CodimOneGridAdapter<typename GridManager::Embeddings>;
    BulkFacetGridAdapter facetGridAdapter(gridManager.getEmbeddings());
    gridGeometry.update(facetGridView, facetGridAdapter);
}

/*!
 * \brief Updates the finite volume grid geometry for the cell-centered schemes.
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

// main program
int main(int argc, char** argv)
{
    using namespace Dumux;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // initialize parameter tree
    Parameters::init(argc, argv);

    // create the grids (from the given grid file)
    using BulkGrid = GetPropType<BulkTypeTag, Properties::Grid>;
    using FacetGrid = GetPropType<FacetTypeTag, Properties::Grid>;
    using GridManager = FacetCouplingGridManager<BulkGrid, FacetGrid>;
    GridManager gridManager;
    gridManager.init();
    gridManager.loadBalance();

    // we compute on the leaf grid views
    const auto& bulkGridView = gridManager.template grid<0>().leafGridView();
    const auto& facetGridView = gridManager.template grid<1>().leafGridView();

    // create the finite volume grid geometries
    using BulkFVGridGeometry = GetPropType<BulkTypeTag, Properties::GridGeometry>;
    using FacetFVGridGeometry = GetPropType<FacetTypeTag, Properties::GridGeometry>;
    auto bulkFvGridGeometry = std::make_shared<BulkFVGridGeometry>(bulkGridView);
    auto facetFvGridGeometry = std::make_shared<FacetFVGridGeometry>(facetGridView);
    updateBulkFVGridGeometry(*bulkFvGridGeometry, gridManager, facetGridView);
    facetFvGridGeometry->update();

    // the coupling mapper
    auto couplingMapper = std::make_shared<typename TestTraits::CouplingMapper>();
    couplingMapper->update(*bulkFvGridGeometry, *facetFvGridGeometry, gridManager.getEmbeddings());

    // the coupling manager
    using CouplingManager = typename TestTraits::CouplingManager;
    auto couplingManager = std::make_shared<CouplingManager>();

    // the problems (boundary conditions)
    using BulkProblem = GetPropType<BulkTypeTag, Properties::Problem>;
    using FacetProblem = GetPropType<FacetTypeTag, Properties::Problem>;
    auto bulkSpatialParams = std::make_shared<typename BulkProblem::SpatialParams>(bulkFvGridGeometry, "Bulk");
    auto bulkProblem = std::make_shared<BulkProblem>(bulkFvGridGeometry, bulkSpatialParams, couplingManager, "Bulk");
    auto facetSpatialParams = std::make_shared<typename FacetProblem::SpatialParams>(facetFvGridGeometry, "Facet");
    auto facetProblem = std::make_shared<FacetProblem>(facetFvGridGeometry, facetSpatialParams, couplingManager, "Facet");

    // the solution vector
    using MDTraits = typename TestTraits::MDTraits;
    using SolutionVector = typename MDTraits::SolutionVector;
    SolutionVector x, xOld;

    static const auto bulkId = typename MDTraits::template SubDomain<0>::Index();
    static const auto facetId = typename MDTraits::template SubDomain<1>::Index();
    x[bulkId].resize(bulkFvGridGeometry->numDofs());
    x[facetId].resize(facetFvGridGeometry->numDofs());
    bulkProblem->applyInitialSolution(x[bulkId]);
    facetProblem->applyInitialSolution(x[facetId]);
    xOld = x;

    // initialize coupling manager
    couplingManager->init(bulkProblem, facetProblem, couplingMapper, x);

    // the grid variables
    using BulkGridVariables = GetPropType<BulkTypeTag, Properties::GridVariables>;
    using FacetGridVariables = GetPropType<FacetTypeTag, Properties::GridVariables>;
    auto bulkGridVariables = std::make_shared<BulkGridVariables>(bulkProblem, bulkFvGridGeometry);
    auto facetGridVariables = std::make_shared<FacetGridVariables>(facetProblem, facetFvGridGeometry);
    bulkGridVariables->init(x[bulkId]);
    facetGridVariables->init(x[facetId]);

    // intialize the vtk output module
    using BulkSolutionVector = std::decay_t<decltype(x[bulkId])>;
    using FacetSolutionVector = std::decay_t<decltype(x[facetId])>;
    VtkOutputModule<BulkGridVariables, BulkSolutionVector> bulkVtkWriter(*bulkGridVariables, x[bulkId], bulkProblem->name(), "Bulk");
    VtkOutputModule<FacetGridVariables, FacetSolutionVector> facetVtkWriter(*facetGridVariables, x[facetId], facetProblem->name(), "Facet");

    // Add model specific output fields
    using BulkIOFields = GetPropType<BulkTypeTag, Properties::IOFields>;
    using FacetIOFields = GetPropType<FacetTypeTag, Properties::IOFields>;
    BulkIOFields::initOutputModule(bulkVtkWriter);
    FacetIOFields::initOutputModule(facetVtkWriter);

    // write initial solution
    bulkVtkWriter.write(0.0);
    facetVtkWriter.write(0.0);

    // instantiate time loop
    auto timeLoop = std::make_shared<TimeLoop<double>>(0.0,
                                                       getParam<double>("TimeLoop.Dt"),
                                                       getParam<double>("TimeLoop.TEnd"));
    // the assembler
    using Assembler = MultiDomainFVAssembler<MDTraits, CouplingManager, DiffMethod::numeric, /*implicit?*/true>;
    auto assembler = std::make_shared<Assembler>( std::make_tuple(bulkProblem, facetProblem),
                                                  std::make_tuple(bulkFvGridGeometry, facetFvGridGeometry),
                                                  std::make_tuple(bulkGridVariables, facetGridVariables),
                                                  couplingManager, timeLoop, xOld);

    // the linear solver
    using LinearSolver = ILU0BiCGSTABBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = Dumux::MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    auto newtonSolver = std::make_shared<NewtonSolver>(assembler, linearSolver, couplingManager);

    // time loop
    timeLoop->start(); do
    {
        // linearize & solve
        newtonSolver->solve(x, *timeLoop);

        // make the new solution the old solution
        bulkGridVariables->advanceTimeStep();
        facetGridVariables->advanceTimeStep();
        xOld = x;

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by the newton solver
        timeLoop->setTimeStepSize(newtonSolver->suggestTimeStepSize(timeLoop->timeStepSize()));

        // write vtk output
        bulkVtkWriter.write(timeLoop->time());
        facetVtkWriter.write(timeLoop->time());
    } while (!timeLoop->finished());

    return 0;
}
