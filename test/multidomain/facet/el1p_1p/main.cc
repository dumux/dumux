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
 * \brief Test for the elastic single-phase model coupled to a
 *        single-phase model in the facet domain.
 */
#include <config.h>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>

#include "problem_bulk_onep.hh"
#include "problem_facet_onep.hh"
#include "problem_bulk_poroelastic.hh"

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
#include <dumux/multidomain/facet/codimonegridadapter.hh>
#include <dumux/multidomain/facet/geomechanics/couplingmanager.hh>
#include <dumux/multidomain/io/vtkoutputmodule.hh>

// obtain/define some types to be used below in the property definitions and in main
class TestTraits
{
    using BulkFlowFVG = Dumux::GetPropType<Dumux::Properties::TTag::OnePBulkTpfa, Dumux::Properties::FVGridGeometry>;
    using FacetFlowFVG = Dumux::GetPropType<Dumux::Properties::TTag::OnePFacetTpfa, Dumux::Properties::FVGridGeometry>;
    using BulkMechFVG = Dumux::GetPropType<Dumux::Properties::TTag::PoroElasticBulk, Dumux::Properties::FVGridGeometry>;
public:
    using MDTraits = Dumux::MultiDomainTraits< Dumux::Properties::TTag::OnePBulkTpfa,
                                               Dumux::Properties::TTag::OnePFacetTpfa,
                                               Dumux::Properties::TTag::PoroElasticBulk >;

    using CouplingMapperFlow = Dumux::FacetCouplingMapper<BulkFlowFVG, FacetFlowFVG>;
    using CouplingMapperMech = Dumux::FacetCouplingMapper<BulkMechFVG, FacetFlowFVG>;
    using CouplingManager = Dumux::FacetCouplingPoroMechanicsCouplingManager<MDTraits, CouplingMapperFlow, CouplingMapperMech>;
};

// set the coupling manager property in the sub-problems
namespace Dumux {
namespace Properties {

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePBulkTpfa> { using type = typename TestTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::OnePFacetTpfa> { using type = typename TestTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::PoroElasticBulk> { using type = typename TestTraits::CouplingManager; };

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
    constexpr auto bulkFlowId = Traits::template SubDomain<0>::Index{};
    constexpr auto facetFlowId = Traits::template SubDomain<1>::Index{};
    constexpr auto bulkMechId = Traits::template SubDomain<2>::Index{};

    // try to create a grid (from the given grid file or the input file)
    using GridManager = FacetCouplingGridManager<Traits::template SubDomain<bulkFlowId>::Grid,
                                                 Traits::template SubDomain<facetFlowId>::Grid>;

    GridManager gridManager;
    gridManager.init();
    gridManager.loadBalance();

    ////////////////////////////////////////////////////////////
    // run stationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid views
    const auto& bulkGridView = gridManager.template grid<0>().leafGridView();
    const auto& facetGridView = gridManager.template grid<1>().leafGridView();

    // create the finite volume grid geometries
    MultiDomainFVGridGeometry<Traits> fvGridGeometry(std::make_tuple(bulkGridView, facetGridView, bulkGridView));
    fvGridGeometry[bulkFlowId].update();
    fvGridGeometry[facetFlowId].update();
    fvGridGeometry[bulkMechId].update(facetGridView, CodimOneGridAdapter<typename GridManager::Embeddings>(gridManager.getEmbeddings()));

    // the coupling manager
    using CouplingManager = typename TestTraits::CouplingManager;
    auto couplingManager = std::make_shared<CouplingManager>();

    // the problems (boundary conditions)
    MultiDomainFVProblem<Traits> problem;

    using BulkFlowProblem = MultiDomainFVProblem<Traits>::template Type<bulkFlowId>;
    using FacetFlowProblem = MultiDomainFVProblem<Traits>::template Type<facetFlowId>;
    using BulkElasticProblem = MultiDomainFVProblem<Traits>::template Type<bulkMechId>;

    auto bulkFlowSpatialParams = std::make_shared<typename BulkFlowProblem::SpatialParams>(fvGridGeometry.get(bulkFlowId), couplingManager, "OnePBulk");
    auto facetFlowSpatialParams = std::make_shared<typename FacetFlowProblem::SpatialParams>(fvGridGeometry.get(facetFlowId), couplingManager, "OnePFacet");
    auto bulkElasticSpatialParams = std::make_shared<typename BulkElasticProblem::SpatialParams>(fvGridGeometry.get(bulkMechId), "ElasticBulk");

    problem.set(std::make_shared<BulkFlowProblem>(fvGridGeometry.get(bulkFlowId), bulkFlowSpatialParams, couplingManager, "OnePBulk"), bulkFlowId);
    problem.set(std::make_shared<FacetFlowProblem>(fvGridGeometry.get(facetFlowId), facetFlowSpatialParams, couplingManager, "OnePFacet"), facetFlowId);
    problem.set(std::make_shared<BulkElasticProblem>(fvGridGeometry.get(bulkMechId), bulkElasticSpatialParams, couplingManager, "ElasticBulk"), bulkMechId);

    // the solution vector
    typename Traits::SolutionVector x;
    problem.applyInitialSolution(x);

    // the coupling mappers
    using CouplingMapperFlow = typename TestTraits::CouplingMapperFlow;
    using CouplingMapperMech = typename TestTraits::CouplingMapperMech;

    auto couplingMapperFlow = std::make_shared<CouplingMapperFlow>();
    auto couplingMapperMech = std::make_shared<CouplingMapperMech>();
    couplingMapperFlow->update(fvGridGeometry[bulkFlowId], fvGridGeometry[facetFlowId], gridManager.getEmbeddings());
    couplingMapperMech->update(fvGridGeometry[bulkMechId], fvGridGeometry[facetFlowId], gridManager.getEmbeddings());

    // initialize the coupling manager
    couplingManager->init(problem.get(bulkFlowId),
                          problem.get(facetFlowId),
                          problem.get(bulkMechId),
                          couplingMapperFlow,
                          couplingMapperMech,
                          x);

    // the grid variables
    using GridVariables = MultiDomainFVGridVariables<Traits>;
    GridVariables gridVars(fvGridGeometry.getTuple(), problem.getTuple());
    gridVars.init(x);

    // intialize the vtk output module
    MultiDomainVtkOutputModule<Traits> vtkWriter;

    using BulkFlowVtkOutputModule = typename MultiDomainVtkOutputModule<Traits>::template Type<bulkFlowId>;
    using FacetFlowVtkOutputModule = typename MultiDomainVtkOutputModule<Traits>::template Type<facetFlowId>;
    using BulkElasticVtkOutputModule = typename MultiDomainVtkOutputModule<Traits>::template Type<bulkMechId>;
    vtkWriter.set(std::make_shared<BulkFlowVtkOutputModule>(gridVars[bulkFlowId], x[bulkFlowId], problem[bulkFlowId].name()), bulkFlowId);
    vtkWriter.set(std::make_shared<FacetFlowVtkOutputModule>(gridVars[facetFlowId], x[facetFlowId], problem[facetFlowId].name()), facetFlowId);
    vtkWriter.set(std::make_shared<BulkElasticVtkOutputModule>(gridVars[bulkMechId], x[bulkMechId], problem[bulkMechId].name(), "ElasticBulk", Dune::VTK::nonconforming), bulkMechId);

    // add additional output
    vtkWriter[bulkFlowId].addField(problem[bulkFlowId].porosities(), "porosity");
    vtkWriter[bulkFlowId].addField(problem[bulkFlowId].permeabilities(), "permeability");
    vtkWriter[facetFlowId].addField(problem[facetFlowId].apertures(), "aperture");
    vtkWriter[facetFlowId].addField(problem[facetFlowId].permeabilities(), "permeability");

    vtkWriter.initDefaultOutputFields();
    vtkWriter.write(0.0);

    // the assembler
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric, /*implicit?*/true>;
    auto assembler = std::make_shared<Assembler>( problem.getTuple(), fvGridGeometry.getTuple(), gridVars.getTuple(), couplingManager);

    // the linear solver
    using LinearSolver = UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = Dumux::MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    auto newtonSolver = std::make_shared<NewtonSolver>(assembler, linearSolver, couplingManager);

    // linearize & solve
    newtonSolver->solve(x);

    // update grid variables for output
    gridVars.update(x);

    // update additional output fields
    problem[bulkFlowId].updateOutputFields(x[bulkFlowId]);
    problem[facetFlowId].updateOutputFields(x[facetFlowId]);

    // write vtk output
    vtkWriter.write(1.0);

    // print parameter usage
    Parameters::print();

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
