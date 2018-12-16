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
#include <dumux/common/defaultusagemessage.hh>

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
#include <dumux/multidomain/io/vtkoutputmodule.hh>

// obtain/define some types to be used below in the property definitions and in main
class TestTraits
{
    using BulkFVG = Dumux::GetPropType<Dumux::Properties::TTag::OnePBulkTpfa, Dumux::Properties::FVGridGeometry>;
    using FacetFVG = Dumux::GetPropType<Dumux::Properties::TTag::OnePFacetTpfa, Dumux::Properties::FVGridGeometry>;
    using EdgeFVG = Dumux::GetPropType<Dumux::Properties::TTag::OnePEdgeTpfa, Dumux::Properties::FVGridGeometry>;
public:
    using MDTraits = Dumux::MultiDomainTraits<Dumux::Properties::TTag::OnePBulkTpfa, Dumux::Properties::TTag::OnePFacetTpfa, Dumux::Properties::TTag::OnePEdgeTpfa>;
    using CouplingMapper = Dumux::FacetCouplingThreeDomainMapper<BulkFVG, FacetFVG, EdgeFVG>;
    using CouplingManager = Dumux::FacetCouplingThreeDomainManager<MDTraits, CouplingMapper>;
};

// set the coupling manager property in the sub-problems
namespace Dumux {
namespace Properties {

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePBulkTpfa> { using type = typename TestTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePFacetTpfa> { using type = typename TestTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePEdgeTpfa> { using type = typename TestTraits::CouplingManager; };

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

    // the multidomain traits and some indices
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
    MultiDomainFVGridGeometry<Traits> fvGridGeometry(std::make_tuple(bulkGridView, facetGridView, edgeGridView));
    fvGridGeometry.update();

    // the coupling manager
    using CouplingManager = typename TestTraits::CouplingManager;
    auto couplingManager = std::make_shared<CouplingManager>();

    // the problems (boundary conditions)
    MultiDomainFVProblem<Traits> problem;

    using BulkProblem = MultiDomainFVProblem<Traits>::template Type<bulkId>;
    using FacetProblem = MultiDomainFVProblem<Traits>::template Type<facetId>;
    using EdgeProblem = MultiDomainFVProblem<Traits>::template Type<edgeId>;

    auto bulkSpatialParams = std::make_shared<typename BulkProblem::SpatialParams>(fvGridGeometry.get(bulkId), "Bulk");
    auto facetSpatialParams = std::make_shared<typename FacetProblem::SpatialParams>(fvGridGeometry.get(facetId), "Facet");
    auto edgeSpatialParams = std::make_shared<typename EdgeProblem::SpatialParams>(fvGridGeometry.get(edgeId), "Edge");

    problem.set(std::make_shared<BulkProblem>(fvGridGeometry.get(bulkId), bulkSpatialParams, couplingManager, "Bulk"), bulkId);
    problem.set(std::make_shared<FacetProblem>(fvGridGeometry.get(facetId), facetSpatialParams, couplingManager, "Facet"), facetId);
    problem.set(std::make_shared<EdgeProblem>(fvGridGeometry.get(edgeId), edgeSpatialParams, couplingManager, "Edge"), edgeId);

    // the solution vector
    typename Traits::SolutionVector x;
    problem.applyInitialSolution(x);

    // the coupling mapper
    using CouplingMapper = typename TestTraits::CouplingMapper;
    auto couplingMapper = std::make_shared<CouplingMapper>();
    couplingMapper->update(fvGridGeometry[bulkId], fvGridGeometry[facetId], fvGridGeometry[edgeId], gridManager.getEmbeddings());

    // initialize the coupling manager
    couplingManager->init(problem.get(bulkId), problem.get(facetId), problem.get(edgeId), couplingMapper, x);

    // the grid variables
    using GridVariables = MultiDomainFVGridVariables<Traits>;
    GridVariables gridVars(fvGridGeometry.getTuple(), problem.getTuple());
    gridVars.init(x);

    // intialize the vtk output module
    const std::array<std::string, 3> vtkOutputNames{{problem[bulkId].name(), problem[facetId].name(), problem[edgeId].name()}};
    MultiDomainVtkOutputModule<Traits> vtkWriter(gridVars.getTuple(), x, vtkOutputNames);
    vtkWriter.initDefaultOutputFields();
    vtkWriter.write(0.0);

    // the assembler
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric, /*implicit?*/true>;
    auto assembler = std::make_shared<Assembler>( problem.getTuple(), fvGridGeometry.getTuple(), gridVars.getTuple(), couplingManager);

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
