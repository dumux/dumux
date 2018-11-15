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

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>

#include "problem_bulk.hh"
#include "problem_facet.hh"
#include "problem_edge.hh"

#include <dumux/assembly/diffmethod.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>

#include <dumux/multidomain/facet/gridmanager.hh>
#include <dumux/multidomain/facet/couplingmapper.hh>
#include <dumux/multidomain/facet/couplingmanager.hh>

#include <dumux/io/vtkoutputmodule.hh>

// obtain/define some types to be used below in the property definitions and in main
class TestTraits
{
    using BulkFVG = Dumux::GetPropType<TTAG(OnePBulkTpfa), Dumux::Properties::FVGridGeometry>;
    using FacetFVG = Dumux::GetPropType<TTAG(OnePFacetTpfa), Dumux::Properties::FVGridGeometry>;
    using EdgeFVG = Dumux::GetPropType<TTAG(OnePEdgeTpfa), Dumux::Properties::FVGridGeometry>;
public:
    using MDTraits = Dumux::MultiDomainTraits<TTAG(OnePBulkTpfa), TTAG(OnePFacetTpfa), TTAG(OnePEdgeTpfa)>;
    using CouplingMapper = Dumux::FacetCouplingThreeDomainMapper<BulkFVG, FacetFVG, EdgeFVG>;
    using CouplingManager = Dumux::FacetCouplingThreeDomainManager<MDTraits, CouplingMapper>;
};

// set the coupling manager property in the sub-problems
namespace Dumux {
namespace Properties {

SET_TYPE_PROP(OnePBulkTpfa, CouplingManager, typename TestTraits::CouplingManager);
SET_TYPE_PROP(OnePFacetTpfa, CouplingManager, typename TestTraits::CouplingManager);
SET_TYPE_PROP(OnePEdgeTpfa, CouplingManager, typename TestTraits::CouplingManager);

} // end namespace Properties
} // end namespace Dumux

int main(int argc, char** argv) try
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

    // try to create a grid (from the given grid file or the input file)
    using BulkProblemTypeTag = TTAG(OnePBulkTpfa);
    using FacetProblemTypeTag = TTAG(OnePFacetTpfa);
    using EdgeProblemTypeTag = TTAG(OnePEdgeTpfa);
    using BulkGrid = GetPropType<BulkProblemTypeTag, Properties::Grid>;
    using FacetGrid = GetPropType<FacetProblemTypeTag, Properties::Grid>;
    using EdgeGrid = GetPropType<EdgeProblemTypeTag, Properties::Grid>;

    using GridManager = FacetCouplingGridManager<BulkGrid, FacetGrid, EdgeGrid>;
    GridManager gridManager;
    gridManager.init();
    gridManager.loadBalance();

    ////////////////////////////////////////////////////////////
    // run stationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid views
    const auto& bulkGridView = gridManager.template grid<0>().leafGridView();
    const auto& facetGridView = gridManager.template grid<1>().leafGridView();
    const auto& edgeGridView = gridManager.template grid<2>().leafGridView();

    // create the finite volume grid geometries
    using BulkFVGridGeometry = GetPropType<BulkProblemTypeTag, Properties::FVGridGeometry>;
    using FacetFVGridGeometry = GetPropType<FacetProblemTypeTag, Properties::FVGridGeometry>;
    using EdgeFVGridGeometry = GetPropType<EdgeProblemTypeTag, Properties::FVGridGeometry>;
    auto bulkFvGridGeometry = std::make_shared<BulkFVGridGeometry>(bulkGridView);
    auto facetFvGridGeometry = std::make_shared<FacetFVGridGeometry>(facetGridView);
    auto edgeFvGridGeometry = std::make_shared<EdgeFVGridGeometry>(edgeGridView);
    bulkFvGridGeometry->update();
    facetFvGridGeometry->update();
    edgeFvGridGeometry->update();

    // the coupling manager
    using CouplingManager = typename TestTraits::CouplingManager;
    auto couplingManager = std::make_shared<CouplingManager>();

    // the problems (boundary conditions)
    using BulkProblem = GetPropType<BulkProblemTypeTag, Properties::Problem>;
    using FacetProblem = GetPropType<FacetProblemTypeTag, Properties::Problem>;
    using EdgeProblem = GetPropType<EdgeProblemTypeTag, Properties::Problem>;
    auto bulkSpatialParams = std::make_shared<typename BulkProblem::SpatialParams>(bulkFvGridGeometry, "Bulk");
    auto bulkProblem = std::make_shared<BulkProblem>(bulkFvGridGeometry, bulkSpatialParams, couplingManager, "Bulk");
    auto facetSpatialParams = std::make_shared<typename FacetProblem::SpatialParams>(facetFvGridGeometry, "Facet");
    auto facetProblem = std::make_shared<FacetProblem>(facetFvGridGeometry, facetSpatialParams, couplingManager, "Facet");
    auto edgeSpatialParams = std::make_shared<typename EdgeProblem::SpatialParams>(edgeFvGridGeometry, "Edge");
    auto edgeProblem = std::make_shared<EdgeProblem>(edgeFvGridGeometry, edgeSpatialParams, couplingManager, "Edge");

    // the solution vector
    using Traits = typename TestTraits::MDTraits;
    using SolutionVector = typename Traits::SolutionVector;
    SolutionVector x;

    static const auto bulkId = Traits::template DomainIdx<0>();
    static const auto facetId = Traits::template DomainIdx<1>();
    static const auto edgeId = Traits::template DomainIdx<2>();
    x[bulkId].resize(bulkFvGridGeometry->numDofs());
    x[facetId].resize(facetFvGridGeometry->numDofs());
    x[edgeId].resize(edgeFvGridGeometry->numDofs());
    bulkProblem->applyInitialSolution(x[bulkId]);
    facetProblem->applyInitialSolution(x[facetId]);
    edgeProblem->applyInitialSolution(x[edgeId]);

    // the coupling mapper
    using CouplingMapper = typename TestTraits::CouplingMapper;
    auto couplingMapper = std::make_shared<CouplingMapper>();
    couplingMapper->update(*bulkFvGridGeometry, *facetFvGridGeometry, *edgeFvGridGeometry, gridManager.getEmbeddings());

    // initialize the coupling manager
    couplingManager->init(bulkProblem, facetProblem, edgeProblem, couplingMapper, x);

    // the grid variables
    using BulkGridVariables = GetPropType<BulkProblemTypeTag, Properties::GridVariables>;
    using FacetGridVariables = GetPropType<FacetProblemTypeTag, Properties::GridVariables>;
    using EdgeGridVariables = GetPropType<EdgeProblemTypeTag, Properties::GridVariables>;
    auto bulkGridVariables = std::make_shared<BulkGridVariables>(bulkProblem, bulkFvGridGeometry);
    auto facetGridVariables = std::make_shared<FacetGridVariables>(facetProblem, facetFvGridGeometry);
    auto edgeGridVariables = std::make_shared<EdgeGridVariables>(edgeProblem, edgeFvGridGeometry);
    bulkGridVariables->init(x[bulkId]);
    facetGridVariables->init(x[facetId]);
    edgeGridVariables->init(x[edgeId]);

    // intialize the vtk output module
    using BulkSolutionVector = std::decay_t<decltype(x[bulkId])>;
    using FacetSolutionVector = std::decay_t<decltype(x[facetId])>;
    using EdgeSolutionVector = std::decay_t<decltype(x[edgeId])>;

    // intialize the vtk output module
    VtkOutputModule<BulkGridVariables, BulkSolutionVector> bulkVtkWriter(*bulkGridVariables, x[bulkId], bulkProblem->name());
    VtkOutputModule<FacetGridVariables, FacetSolutionVector> facetVtkWriter(*facetGridVariables, x[facetId],  facetProblem->name());
    VtkOutputModule<EdgeGridVariables, EdgeSolutionVector> edgeVtkWriter(*edgeGridVariables, x[edgeId], edgeProblem->name());

    // Add model specific output fields
    using BulkIOFields = GetPropType<BulkProblemTypeTag, Properties::IOFields>;
    using FacetIOFields = GetPropType<FacetProblemTypeTag, Properties::IOFields>;
    using EdgeIOFields = GetPropType<EdgeProblemTypeTag, Properties::IOFields>;
    BulkIOFields::initOutputModule(bulkVtkWriter);
    FacetIOFields::initOutputModule(facetVtkWriter);
    EdgeIOFields::initOutputModule(edgeVtkWriter);
    bulkVtkWriter.write(0.0);
    facetVtkWriter.write(0.0);
    edgeVtkWriter.write(0.0);

    // the assembler
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric, /*implicit?*/true>;
    auto assembler = std::make_shared<Assembler>( std::make_tuple(bulkProblem, facetProblem, edgeProblem),
                                                  std::make_tuple(bulkFvGridGeometry, facetFvGridGeometry, edgeFvGridGeometry),
                                                  std::make_tuple(bulkGridVariables, facetGridVariables, edgeGridVariables),
                                                  couplingManager);

    // the linear solver
    using LinearSolver = ILU0BiCGSTABBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = Dumux::MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    auto newtonSolver = std::make_shared<NewtonSolver>(assembler, linearSolver, couplingManager);

    // linearize & solve
    newtonSolver->solve(x);

    // update grid variables for output
    bulkGridVariables->update(x[bulkId]);
    facetGridVariables->update(x[facetId]);
    edgeGridVariables->update(x[edgeId]);

    // write vtk output
    bulkVtkWriter.write(1.0);
    facetVtkWriter.write(1.0);
    edgeVtkWriter.write(1.0);

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
