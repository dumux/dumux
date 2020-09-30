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

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include "problem_bulk.hh"
#include "problem_facet.hh"
#include "problem_edge.hh"

#include <dumux/assembly/diffmethod.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/fvgridgeometry.hh>
#include <dumux/multidomain/fvproblem.hh>
#include <dumux/multidomain/fvgridvariables.hh>

#include <dumux/multidomain/facet/gridmanager.hh>
#include <dumux/multidomain/facet/couplingmapper.hh>
#include <dumux/multidomain/facet/couplingmanager.hh>
#include <dumux/multidomain/facet/codimonegridadapter.hh>
#include <dumux/multidomain/io/vtkoutputmodule.hh>

// obtain/define some types to be used below in the property definitions and in main
template<class BulkTypeTag, class FacetTypeTag, class EdgeTypeTag>
class TestTraits
{
    using BulkFVG = Dumux::GetPropType<BulkTypeTag, Dumux::Properties::GridGeometry>;
    using FacetFVG = Dumux::GetPropType<FacetTypeTag, Dumux::Properties::GridGeometry>;
    using EdgeFVG = Dumux::GetPropType<EdgeTypeTag, Dumux::Properties::GridGeometry>;
public:
    using MDTraits = Dumux::MultiDomainTraits<BulkTypeTag, FacetTypeTag, EdgeTypeTag>;
    using CouplingMapper = Dumux::FacetCouplingThreeDomainMapper<BulkFVG, FacetFVG, EdgeFVG>;
    using CouplingManager = Dumux::FacetCouplingThreeDomainManager<MDTraits, CouplingMapper>;
};

// set the coupling manager property in the sub-problems
namespace Dumux {
namespace Properties {

using TpfaTraits = TestTraits<TTag::OnePBulkTpfa, TTag::OnePFacetTpfa, TTag::OnePEdgeTpfa>;
using MpfaTraits = TestTraits<TTag::OnePBulkMpfa, TTag::OnePFacetMpfa, TTag::OnePEdgeMpfa>;
using BoxTraits = TestTraits<TTag::OnePBulkBox, TTag::OnePFacetBox, TTag::OnePEdgeBox>;

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePBulkTpfa> { using type = typename TpfaTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePFacetTpfa> { using type = typename TpfaTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePEdgeTpfa> { using type = typename TpfaTraits::CouplingManager; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePBulkMpfa> { using type = typename MpfaTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePFacetMpfa> { using type = typename MpfaTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePEdgeMpfa> { using type = typename MpfaTraits::CouplingManager; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePBulkBox> { using type = typename BoxTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePFacetBox> { using type = typename BoxTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePEdgeBox> { using type = typename BoxTraits::CouplingManager; };

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
          class LowDimGridView,
          std::enable_if_t<GridGeometry::discMethod == Dumux::DiscretizationMethod::box, int> = 0 >
void updateFVGridGeometry(GridGeometry& gridGeometry,
                          const GridManager& gridManager,
                          const LowDimGridView& lowDimGridView)
{
    static constexpr int higherGridId = int(GridGeometry::GridView::dimension) == 3 ? 0 : 1;
    using BulkFacetGridAdapter = Dumux::CodimOneGridAdapter<typename GridManager::Embeddings, higherGridId, higherGridId+1>;
    BulkFacetGridAdapter facetGridAdapter(gridManager.getEmbeddings());
    gridGeometry.update(lowDimGridView, facetGridAdapter);
}

/*!
 * \brief Updates the finite volume grid geometry for the cell-centered schemes.
 */
template< class GridGeometry,
          class GridManager,
          class LowDimGridView,
          std::enable_if_t<GridGeometry::discMethod != Dumux::DiscretizationMethod::box, int> = 0 >
void updateFVGridGeometry(GridGeometry& gridGeometry,
                          const GridManager& gridManager,
                          const LowDimGridView& lowDimGridView)
{
    gridGeometry.update();
}

int main(int argc, char** argv)
{
    using namespace Dumux;

    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

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
    using ThisTestTraits = TestTraits<BulkTypeTag, FacetTypeTag, EdgeTypeTag>;
    using Traits = typename ThisTestTraits::MDTraits;
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
    MultiDomainFVGridGeometry<Traits> gridGeometry(std::make_tuple(bulkGridView, facetGridView, edgeGridView));
    updateFVGridGeometry(gridGeometry[bulkId], gridManager, facetGridView);
    updateFVGridGeometry(gridGeometry[facetId], gridManager, edgeGridView);
    gridGeometry[edgeId].update();

    // the coupling manager
    using CouplingManager = typename ThisTestTraits::CouplingManager;
    auto couplingManager = std::make_shared<CouplingManager>();

    // the problems (boundary conditions)
    MultiDomainFVProblem<Traits> problem;

    using BulkProblem = MultiDomainFVProblem<Traits>::template Type<bulkId>;
    using FacetProblem = MultiDomainFVProblem<Traits>::template Type<facetId>;
    using EdgeProblem = MultiDomainFVProblem<Traits>::template Type<edgeId>;

    auto bulkSpatialParams = std::make_shared<typename BulkProblem::SpatialParams>(gridGeometry.get(bulkId), "Bulk");
    auto facetSpatialParams = std::make_shared<typename FacetProblem::SpatialParams>(gridGeometry.get(facetId), "Facet");
    auto edgeSpatialParams = std::make_shared<typename EdgeProblem::SpatialParams>(gridGeometry.get(edgeId), "Edge");

    problem.set(std::make_shared<BulkProblem>(gridGeometry.get(bulkId), bulkSpatialParams, couplingManager, "Bulk"), bulkId);
    problem.set(std::make_shared<FacetProblem>(gridGeometry.get(facetId), facetSpatialParams, couplingManager, "Facet"), facetId);
    problem.set(std::make_shared<EdgeProblem>(gridGeometry.get(edgeId), edgeSpatialParams, couplingManager, "Edge"), edgeId);

    // the solution vector
    typename Traits::SolutionVector x;
    problem.applyInitialSolution(x);

    // the coupling mapper
    using CouplingMapper = typename ThisTestTraits::CouplingMapper;
    auto couplingMapper = std::make_shared<CouplingMapper>();
    couplingMapper->update(gridGeometry[bulkId], gridGeometry[facetId], gridGeometry[edgeId], gridManager.getEmbeddings());

    // initialize the coupling manager
    couplingManager->init(problem.get(bulkId), problem.get(facetId), problem.get(edgeId), couplingMapper, x);

    // the grid variables
    using GridVariables = MultiDomainFVGridVariables<Traits>;
    GridVariables gridVars(gridGeometry.getTuple(), problem.getTuple());
    gridVars.init(x);

    // intialize the vtk output module
    const std::array<std::string, 3> vtkOutputNames{{problem[bulkId].name(), problem[facetId].name(), problem[edgeId].name()}};
    MultiDomainVtkOutputModule<Traits> vtkWriter(gridVars.getTuple(), x, vtkOutputNames);
    vtkWriter.initDefaultOutputFields();
    vtkWriter.write(0.0);

    // the assembler
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric, /*implicit?*/true>;
    auto assembler = std::make_shared<Assembler>( problem.getTuple(), gridGeometry.getTuple(), gridVars.getTuple(), couplingManager);

    // the linear solver
    using LinearSolver = ILU0BiCGSTABBackend;
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
