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

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>

#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/fvgridgeometry.hh>
#include <dumux/multidomain/fvproblem.hh>
#include <dumux/multidomain/fvgridvariables.hh>

#include <dumux/multidomain/facet/gridmanager.hh>
#include <dumux/multidomain/facet/codimonegridadapter.hh>
#include <dumux/multidomain/io/vtkoutputmodule.hh>

#include "properties.hh"

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
          class LowDimGridView,
          std::enable_if_t<GridGeometry::discMethod == Dumux::DiscretizationMethods::box, int> = 0 >
auto gridGeometryArgs(const typename GridGeometry::GridView& gridView,
                      const GridManager& gridManager,
                      const LowDimGridView& lowDimGridView)
{
    static constexpr int higherGridId = int(GridGeometry::GridView::dimension) == 3 ? 0 : 1;
    using BulkFacetGridAdapter = Dumux::CodimOneGridAdapter<typename GridManager::Embeddings, higherGridId, higherGridId+1>;
    BulkFacetGridAdapter facetGridAdapter(gridManager.getEmbeddings());
    return std::make_tuple(gridView, lowDimGridView, std::move(facetGridAdapter));
}

/*!
 * \brief Updates the finite volume grid geometry for the cell-centered schemes.
 */
template< class GridGeometry,
          class GridManager,
          class LowDimGridView,
          std::enable_if_t<GridGeometry::discMethod != Dumux::DiscretizationMethods::box, int> = 0 >
const typename GridGeometry::GridView&
gridGeometryArgs(const typename GridGeometry::GridView& gridView,
                 const GridManager& gridManager,
                 const LowDimGridView& lowDimGridView)
{ return gridView; }


int main(int argc, char** argv)
{
    using namespace Dumux;

    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // initialize parameter tree
    Parameters::init(argc, argv);

    // The sub-problem type tags
    using BulkTypeTag = BULKTYPETAG;
    using FacetTypeTag = FACETTYPETAG;
    using EdgeTypeTag = EDGETYPETAG;

    // the multidomain traits and some indices
    using TestTraits = Properties::TestTraits<BulkTypeTag, FacetTypeTag, EdgeTypeTag>;
    using Traits = typename TestTraits::MDTraits;
    constexpr auto bulkId = Traits::template SubDomain<0>::Index{};
    constexpr auto facetId = Traits::template SubDomain<1>::Index{};
    constexpr auto edgeId = Traits::template SubDomain<2>::Index{};

    // try to create a grid (from the given grid file or the input file)
    FacetCouplingGridManager<Traits::template SubDomain<bulkId>::Grid,
                             Traits::template SubDomain<facetId>::Grid,
                             Traits::template SubDomain<edgeId>::Grid> gridManager;

    gridManager.init();
    gridManager.loadBalance();

    ////////////////////////////////////////////////////////////
    // run stationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid views
    const auto& bulkGridView = gridManager.template grid<bulkId>().leafGridView();
    const auto& facetGridView = gridManager.template grid<facetId>().leafGridView();
    const auto& edgeGridView = gridManager.template grid<edgeId>().leafGridView();

    // create the finite volume grid geometries
    using BulkGridGeometry = typename MultiDomainFVGridGeometry<Traits>::template Type<bulkId>;
    using FacetGridGeometry = typename MultiDomainFVGridGeometry<Traits>::template Type<facetId>;

    MultiDomainFVGridGeometry<Traits> gridGeometry(
        gridGeometryArgs<BulkGridGeometry>(bulkGridView, gridManager, facetGridView),
        gridGeometryArgs<FacetGridGeometry>(facetGridView, gridManager, edgeGridView),
        edgeGridView
    );

    // the coupling manager
    using CouplingManager = typename TestTraits::CouplingManager;
    auto couplingManager = std::make_shared<CouplingManager>();

    // the problems (boundary conditions)
    using BulkProblem = MultiDomainFVProblem<Traits>::template Type<bulkId>;
    using FacetProblem = MultiDomainFVProblem<Traits>::template Type<facetId>;
    using EdgeProblem = MultiDomainFVProblem<Traits>::template Type<edgeId>;

    auto bulkSpatialParams = std::make_shared<typename BulkProblem::SpatialParams>(gridGeometry.get(bulkId), "Bulk");
    auto facetSpatialParams = std::make_shared<typename FacetProblem::SpatialParams>(gridGeometry.get(facetId), "Facet");
    auto edgeSpatialParams = std::make_shared<typename EdgeProblem::SpatialParams>(gridGeometry.get(edgeId), "Edge");

    auto bulkProblem = std::make_shared<BulkProblem>(gridGeometry.get(bulkId), bulkSpatialParams, couplingManager, "Bulk");
    auto facetProblem = std::make_shared<FacetProblem>(gridGeometry.get(facetId), facetSpatialParams, couplingManager, "Facet");
    auto edgeProblem = std::make_shared<EdgeProblem>(gridGeometry.get(edgeId), edgeSpatialParams, couplingManager, "Edge");

    MultiDomainFVProblem<Traits> problem(std::make_tuple(bulkProblem, facetProblem, edgeProblem));

    // the solution vector
    typename Traits::SolutionVector x;
    problem.applyInitialSolution(x);

    // the coupling mapper
    using CouplingMapper = typename TestTraits::CouplingMapper;
    auto couplingMapper = std::make_shared<CouplingMapper>();
    couplingMapper->update(gridGeometry[bulkId], gridGeometry[facetId], gridGeometry[edgeId], gridManager.getEmbeddings());

    // initialize the coupling manager
    couplingManager->init(problem.get(bulkId), problem.get(facetId), problem.get(edgeId), couplingMapper, x);

    // the grid variables
    using GridVariables = MultiDomainFVGridVariables<Traits>;
    GridVariables gridVars(gridGeometry.asTuple(), problem.asTuple());
    gridVars.init(x);

    // initialize the vtk output module
    const std::array<std::string, 3> vtkOutputNames{{problem[bulkId].name(), problem[facetId].name(), problem[edgeId].name()}};
    MultiDomainVtkOutputModule<Traits> vtkWriter(gridVars.asTuple(), x, vtkOutputNames);
    vtkWriter.initDefaultOutputFields();
    vtkWriter.write(0.0);

    // the assembler
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric, /*implicit?*/true>;
    auto assembler = std::make_shared<Assembler>( problem.asTuple(), gridGeometry.asTuple(), gridVars.asTuple(), couplingManager);

    // the linear solver
    using LinearSolver = ILUBiCGSTABIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = Dumux::MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    auto newtonSolver = std::make_shared<NewtonSolver>(assembler, linearSolver, couplingManager);

    // linearize & solve
    newtonSolver->solve(x);

    // update grid variables for output
    gridVars.update(x);

    // write vtk output
    vtkWriter.write(1.0);

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/false);

    return 0;

}
