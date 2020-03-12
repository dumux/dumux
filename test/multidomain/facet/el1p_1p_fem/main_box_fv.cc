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
 * \brief Test for the elastic single-phase model coupled to a
 *        single-phase model in the facet domain, together with
 *        a mechanical model for the deformations and a contact
 *        problem on the fracture facets.
 * \note This test uses the box scheme for the mechanical deformations and
 *       a finite volume approach for the flow field.
 */
#include <config.h>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>

#include "problem_bulk_onep.hh"
#include "problem_facet_onep.hh"
#include "problem_bulk_poroelastic.hh"
#include "problem_lagrangemp.hh"

#include <dumux/io/grid/gridmanager.hh>
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
#include <dumux/multidomain/facet/enrichedgridmanager.hh>
#include <dumux/multidomain/facet/vertexmapper.hh>
#include <dumux/multidomain/facet/geomechanics/couplingmanager.hh>
#include <dumux/multidomain/io/vtkoutputmodule.hh>

// obtain the flow type tags from CMakeLists.txt
using BulkFlowTypeTag = Dumux::Properties::TTag::BULKFLOWTYPETAG;
using FacetFlowTypeTag = Dumux::Properties::TTag::FACETFLOWTYPETAG;

// obtain/define some types to be used below in the property definitions and in main
class TestTraits
{
    using BulkFlowFVG = Dumux::GetPropType<BulkFlowTypeTag, Dumux::Properties::FVGridGeometry>;
    using FacetFlowFVG = Dumux::GetPropType<FacetFlowTypeTag, Dumux::Properties::FVGridGeometry>;
    using BulkMechFVG = Dumux::GetPropType<Dumux::Properties::TTag::PoroElasticBulkBox, Dumux::Properties::FVGridGeometry>;
public:
    using MDTraits = Dumux::MultiDomainTraits< BulkFlowTypeTag,
                                               FacetFlowTypeTag,
                                               Dumux::Properties::TTag::PoroElasticBulkBox,
                                               Dumux::Properties::TTag::LagrangeFacet>;

    using CouplingMapperFlow = Dumux::FacetCouplingMapper<BulkFlowFVG, FacetFlowFVG>;
    using CouplingMapperMech = Dumux::FacetCouplingMapper<BulkMechFVG, FacetFlowFVG>;
    using CouplingManager = Dumux::FacetCouplingPoroMechanicsCouplingManager<MDTraits, CouplingMapperFlow, CouplingMapperMech>;
};

// set the coupling manager property in the sub-problems
namespace Dumux {
namespace Properties {

template<class TypeTag>
struct CouplingManager<TypeTag, BulkFlowTypeTag> { using type = typename TestTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, FacetFlowTypeTag> { using type = typename TestTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::PoroElasticBulkBox> { using type = typename TestTraits::CouplingManager; };
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::LagrangeFacet> { using type = typename TestTraits::CouplingManager; };

} // end namespace Properties
} // end namespace Dumux

//! brief Updates the finite volume grid geometry for the box-facet coupling scheme.
template< class FVGridGeometry, class GridAdapter, class FacetGridView,
          std::enable_if_t<FVGridGeometry::discMethod == Dumux::DiscretizationMethod::box, int> = 0 >
void updateBulkFlowFVGridGeometry(FVGridGeometry& fvGridGeometry,
                                  const FacetGridView& facetGridView,
                                  const GridAdapter& facetGridAdapter)
{ fvGridGeometry.update(facetGridView, facetGridAdapter); }

//! Updates the finite volume grid geometry for the cell-centered facet coupling schemes.
template< class FVGridGeometry, class GridAdapter, class FacetGridView,
          std::enable_if_t<FVGridGeometry::discMethod != Dumux::DiscretizationMethod::box, int> = 0 >
void updateBulkFlowFVGridGeometry(FVGridGeometry& fvGridGeometry,
                                  const FacetGridView& facetGridView,
                                  const GridAdapter& facetGridAdapter)
{ fvGridGeometry.update(); }

///////////////////
// main function //
///////////////////
int main(int argc, char** argv) try
{
    using namespace Dumux;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // initialize parameter tree
    Parameters::init(argc, argv);

    // the multidomain traits and domain indices
    using Traits = typename TestTraits::MDTraits;
    using CouplingManager = typename TestTraits::CouplingManager;
    constexpr auto bulkFlowId = CouplingManager::matrixFlowId;
    constexpr auto facetFlowId = CouplingManager::facetFlowId;
    constexpr auto bulkMechId = CouplingManager::mechanicsId;
    constexpr auto lagrangeId = CouplingManager::lagrangeId;

    //////////////////////
    // Create the grids //
    //////////////////////
    using BulkFlowGrid = Traits::template SubDomain<bulkFlowId>::Grid;
    using FacetFlowGrid = Traits::template SubDomain<facetFlowId>::Grid;
    using GridManager = FacetCouplingGridManager<BulkFlowGrid, FacetFlowGrid>;
    GridManager gridManager;
    gridManager.init();
    gridManager.loadBalance();

    const auto& bulkGridView = gridManager.template grid<0>().leafGridView();
    const auto& facetGridView = gridManager.template grid<1>().leafGridView();

    //////////////////////////////////////////////////////
    // Create the grid used for the lagrange multiplier //
    //////////////////////////////////////////////////////
    using LagrangeGrid = Traits::template SubDomain<lagrangeId>::Grid;
    Dumux::GridManager<LagrangeGrid> lagrangeGridManager;
    lagrangeGridManager.init("Lagrange");
    lagrangeGridManager.loadBalance();

    const auto& lagrangeGridView = lagrangeGridManager.grid().leafGridView();

    // Create the grid geometries of all domains
    using MDGridGeometry = MultiDomainFVGridGeometry<Traits>;
    MDGridGeometry fvGridGeometry;

    // standard construction for bulkFlow and facet flow
    fvGridGeometry.set(std::make_shared<typename MDGridGeometry::template Type<bulkFlowId>>(bulkGridView), bulkFlowId);
    fvGridGeometry.set(std::make_shared<typename MDGridGeometry::template Type<facetFlowId>>(facetGridView), facetFlowId);
    fvGridGeometry.set(std::make_shared<typename MDGridGeometry::template Type<bulkMechId>>(bulkGridView), bulkMechId);

    // FEM-type construction for lagrange multiplier sub-domain
    using LagrangeBasis = typename MDGridGeometry::template Type<lagrangeId>::AnsatzSpaceBasis;
    auto lagrangeBasis = std::make_shared<LagrangeBasis>(lagrangeGridView);
    fvGridGeometry.set(std::make_shared<typename MDGridGeometry::template Type<lagrangeId>>(lagrangeBasis), lagrangeId);

    // update the grid geometries
    using Embeddings = typename GridManager::Embeddings;
    using GridAdapter = CodimOneGridAdapter<Embeddings>;
    GridAdapter bulkFacetGridAdapter(gridManager.getEmbeddings());

    fvGridGeometry[bulkMechId].update(facetGridView, bulkFacetGridAdapter);
    fvGridGeometry[facetFlowId].update();
    fvGridGeometry[lagrangeId].update();

    // update bulk flow grid geometry depending on the scheme used (box requires special update)
    updateBulkFlowFVGridGeometry(fvGridGeometry[bulkFlowId], facetGridView, bulkFacetGridAdapter);

    // the coupling manager
    using CouplingManager = typename TestTraits::CouplingManager;
    auto couplingManager = std::make_shared<CouplingManager>();

    // the problems (boundary conditions)
    MultiDomainFVProblem<Traits> problem;

    using BulkFlowProblem = MultiDomainFVProblem<Traits>::template Type<bulkFlowId>;
    using FacetFlowProblem = MultiDomainFVProblem<Traits>::template Type<facetFlowId>;
    using BulkElasticProblem = MultiDomainFVProblem<Traits>::template Type<bulkMechId>;
    using LagrangeProblem = MultiDomainFVProblem<Traits>::template Type<lagrangeId>;

    auto bulkFlowSpatialParams = std::make_shared<typename BulkFlowProblem::SpatialParams>(fvGridGeometry.get(bulkFlowId), couplingManager, "OnePBulk");
    auto facetFlowSpatialParams = std::make_shared<typename FacetFlowProblem::SpatialParams>(fvGridGeometry.get(facetFlowId), couplingManager, "OnePFacet");
    auto bulkElasticSpatialParams = std::make_shared<typename BulkElasticProblem::SpatialParams>(fvGridGeometry.get(bulkMechId), "ElasticBulk");

    problem.set(std::make_shared<BulkFlowProblem>(fvGridGeometry.get(bulkFlowId), bulkFlowSpatialParams, couplingManager, "OnePBulk"), bulkFlowId);
    problem.set(std::make_shared<FacetFlowProblem>(fvGridGeometry.get(facetFlowId), facetFlowSpatialParams, couplingManager, "OnePFacet"), facetFlowId);
    problem.set(std::make_shared<BulkElasticProblem>(fvGridGeometry.get(bulkMechId), bulkElasticSpatialParams, couplingManager, "ElasticBulk"), bulkMechId);
    problem.set(std::make_shared<LagrangeProblem>(fvGridGeometry.get(lagrangeId), couplingManager, "Lagrange"), lagrangeId);

    // the coupling mappers
    using CouplingMapperFlow = typename TestTraits::CouplingMapperFlow;
    using CouplingMapperMech = typename TestTraits::CouplingMapperMech;

    auto couplingMapperFlow = std::make_shared<CouplingMapperFlow>();
    auto couplingMapperMech = std::make_shared<CouplingMapperMech>();

    couplingMapperFlow->update(fvGridGeometry[bulkFlowId], fvGridGeometry[facetFlowId], gridManager.getEmbeddings());
    couplingMapperMech->update(fvGridGeometry[bulkMechId], fvGridGeometry[facetFlowId], gridManager.getEmbeddings());

    // the solution vector
    typename Traits::SolutionVector x;

    // initialize the coupling manager
    couplingManager->init(problem.get(bulkFlowId), problem.get(facetFlowId), problem.get(bulkMechId), problem.get(lagrangeId),
                          couplingMapperFlow, couplingMapperMech, x);

    // initial values
    problem.applyInitialSolution(x);
    auto xOld = x;

    // the grid variables
    using GridVariables = MultiDomainFVGridVariables<Traits>;
    GridVariables gridVars(fvGridGeometry.getTuple(), problem.getTuple());
    gridVars.init(x);

    // intialize the vtk output modules
    using BulkFlowGridVariables = typename GridVariables::template Type<bulkFlowId>;
    using FacetFlowGridVariables = typename GridVariables::template Type<facetFlowId>;

    using BulkFlowVtkOutputModule = VtkOutputModule<BulkFlowGridVariables, typename Traits::template SubDomain<bulkFlowId>::SolutionVector>;
    using FacetFlowVtkOutputModule = VtkOutputModule<FacetFlowGridVariables, typename Traits::template SubDomain<facetFlowId>::SolutionVector>;

    static constexpr bool bulkFlowIsBox = MDGridGeometry::template Type<bulkFlowId>::discMethod == DiscretizationMethod::box;
    const auto bulkOutputType = bulkFlowIsBox ? Dune::VTK::nonconforming : Dune::VTK::conforming;
    BulkFlowVtkOutputModule bulkFlowVtkWriter(gridVars[bulkFlowId], x[bulkFlowId], problem[bulkFlowId].name(), "", bulkOutputType);
    FacetFlowVtkOutputModule facetFlowVtkWriter(gridVars[facetFlowId], x[facetFlowId], problem[facetFlowId].name());

    // add additional output
    facetFlowVtkWriter.addField(problem[facetFlowId].apertures(), "aperture", FacetFlowVtkOutputModule::FieldType::element);
    facetFlowVtkWriter.addField(problem[facetFlowId].permeabilities(), "permeability", FacetFlowVtkOutputModule::FieldType::element);

    using BulkFlowIOFields = GetPropType<typename Traits::template SubDomain<bulkFlowId>::TypeTag, Properties::IOFields>;
    using FacetFlowIOFields = GetPropType<typename Traits::template SubDomain<facetFlowId>::TypeTag, Properties::IOFields>;
    BulkFlowIOFields::initOutputModule(bulkFlowVtkWriter);
    FacetFlowIOFields::initOutputModule(facetFlowVtkWriter);

    ////////////////////////////////////
    // VTK output for lagrange domain //
    ////////////////////////////////////
    using Vector = Dune::FieldVector<double, LagrangeGrid::dimensionworld>;
    using LagrangeGridView = typename LagrangeGrid::LeafGridView;
    using LagrangeElementMapper = typename MDGridGeometry::template Type<lagrangeId>::ElementMapper;
    using P0LagrangeFunction = Vtk::VectorP0VTKFunction<LagrangeGridView, LagrangeElementMapper, std::vector<Vector>>;

    auto lagrangeWriter = std::make_shared< Dune::VTKWriter<LagrangeGridView> >(lagrangeGridView);
    Dune::VTKSequenceWriter<LagrangeGridView> lagrangeSequenceWriter(lagrangeWriter, getParamFromGroup<std::string>("Lagrange", "Problem.Name"));

    const double initialGap = getParam<double>("SpatialParams.InitialGap");
    std::vector<double> aperture(lagrangeGridView.size(0), initialGap);
    std::vector<double> deltaUt(lagrangeGridView.size(1));
    std::vector<double> deltaUn(lagrangeGridView.size(1));
    std::vector<double> fricCoeff(lagrangeGridView.size(0));
    std::vector<Vector> deltaU(lagrangeGridView.size(0));
    std::vector<Vector> normalTraction(lagrangeGridView.size(0));
    std::vector<Vector> tangentialTraction(lagrangeGridView.size(0));

    auto contactTractions = x[lagrangeId];
    auto gfContactTraction = Dune::Functions::makeDiscreteGlobalBasisFunction<Vector>(*lagrangeBasis, contactTractions);
    auto gfNormTraction = std::make_shared<P0LagrangeFunction>(lagrangeGridView, fvGridGeometry[lagrangeId].elementMapper(), normalTraction, "normalTraction", 2);
    auto gfTangTraction = std::make_shared<P0LagrangeFunction>(lagrangeGridView, fvGridGeometry[lagrangeId].elementMapper(), tangentialTraction, "tangentialTraction", 2);

    lagrangeWriter->addCellData(gfContactTraction, Dune::VTK::FieldInfo("contactTraction", Dune::VTK::FieldInfo::Type::vector, 2));
    lagrangeWriter->addCellData(gfNormTraction);
    lagrangeWriter->addCellData(gfTangTraction);

    auto deltaUFunction = std::make_shared<P0LagrangeFunction>(lagrangeGridView, fvGridGeometry[lagrangeId].elementMapper(), deltaU, "deltaU", 2);
    lagrangeWriter->addCellData(deltaUFunction);
    lagrangeWriter->addVertexData(deltaUt, "delta_Ut");
    lagrangeWriter->addVertexData(deltaUn, "delta_Un");
    lagrangeWriter->addCellData(fricCoeff, "fricCoeff");
    lagrangeWriter->addCellData(aperture, "gap");

    //////////////////////////////////////
    // VTK output for mechanical domain //
    //////////////////////////////////////
    using MechanicalGridView = typename BulkFlowGrid::LeafGridView;
    using MechElementMapper = typename MDGridGeometry::template Type<bulkMechId>::ElementMapper;
    using MechVertexMapper = typename MDGridGeometry::template Type<bulkMechId>::VertexMapper;
    using P0MechanicsFunction = Vtk::VectorP0VTKFunction<MechanicalGridView, MechElementMapper, std::vector<Vector>>;
    using P1MechanicsFunction = Vtk::VectorP1VTKFunction<MechanicalGridView, MechVertexMapper, std::decay_t<decltype(x[bulkMechId])>>;

    auto mechWriter = std::make_shared<Dune::VTKWriter<MechanicalGridView>>(bulkGridView, Dune::VTK::nonconforming);
    Dune::VTKSequenceWriter<MechanicalGridView> mechSequenceWriter(mechWriter, getParamFromGroup<std::string>("ElasticBulk", "Problem.Name"));

    auto uFunction = std::make_shared<P1MechanicsFunction>(bulkGridView, fvGridGeometry[bulkMechId].vertexMapper(), x[bulkMechId], "u", 2);
    auto sigmaXFunction = std::make_shared<P0MechanicsFunction>(bulkGridView, fvGridGeometry[bulkMechId].elementMapper(), problem[bulkMechId].sigma_x(), "Sigma_x", 2);
    auto sigmaYFunction = std::make_shared<P0MechanicsFunction>(bulkGridView, fvGridGeometry[bulkMechId].elementMapper(), problem[bulkMechId].sigma_y(), "Sigma_y", 2);

    mechWriter->addVertexData(uFunction);
    mechWriter->addCellData(sigmaXFunction);
    mechWriter->addCellData(sigmaYFunction);

    // the assembler
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric, /*implicit?*/true>;
    auto assembler = std::make_shared<Assembler>( problem.getTuple(), fvGridGeometry.getTuple(), gridVars.getTuple(), couplingManager);

    // the linear solver
    using LinearSolver = UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = Dumux::MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    auto newtonSolver = std::make_shared<NewtonSolver>(assembler, linearSolver, couplingManager);

    newtonSolver->solve(x);
    gridVars.update(x);

    // update lagrange output fields
    contactTractions = x[lagrangeId];
    for (const auto& element : elements(lagrangeGridView))
    {
        const auto& eg = element.geometry();
        const auto eIdx = fvGridGeometry[lagrangeId].elementMapper().index(element);
        const auto vIdx1 = fvGridGeometry[lagrangeId].vertexMapper().subIndex(element, 0, 1);
        const auto vIdx2 = fvGridGeometry[lagrangeId].vertexMapper().subIndex(element, 1, 1);

        // flip sign if contact surface is not defined in "positive" coord direction
        const auto& contactSurface = couplingManager->getContactSurfaceSegment(element);
        const auto& normal = contactSurface.getBasisVector(1);
        const auto& traction = x[lagrangeId][eIdx];

        normalTraction[eIdx] = normal;
        normalTraction[eIdx] *= traction*normal;

        tangentialTraction[eIdx] = normalTraction[eIdx];
        tangentialTraction[eIdx] *= -1.0;
        tangentialTraction[eIdx] += traction;

        aperture[eIdx] = couplingManager->computeAperture(element, eg.center(), initialGap);
        deltaU[eIdx] = couplingManager->computeDisplacementJump(element, eg.center());
        fricCoeff[eIdx] = problem[lagrangeId].frictionCoefficient(element, eg.center());

        deltaUt[vIdx1] = couplingManager->computeTangentialDisplacementJump(element, eg.corner(0)).two_norm();
        deltaUt[vIdx2] = couplingManager->computeTangentialDisplacementJump(element, eg.corner(1)).two_norm();
        deltaUn[vIdx1] = couplingManager->computeNormalDisplacementJump(element, eg.corner(0));
        deltaUn[vIdx2] = couplingManager->computeNormalDisplacementJump(element, eg.corner(1));

        // define unique orientation
        // contact surfaces might be defined with
        // different orientation for the segments
        if ( (normal*Vector({1.0, 0.0})) < 0.0)
        {
            normalTraction[eIdx] *= -1.0;
            tangentialTraction[eIdx] *= -1.0;
            contactTractions[eIdx] *= -1.0;
            deltaU[eIdx] *= -1.0;
        }
    }

    // update stresses
    problem[bulkMechId].updateOutputFields(gridVars[bulkMechId], *assembler, x[bulkMechId], bulkMechId);

    // write vtk output
    bulkFlowVtkWriter.write(1.0);
    facetFlowVtkWriter.write(1.0);
    lagrangeSequenceWriter.write(1.0);
    mechSequenceWriter.write(1.0);

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
catch (std::exception& e)
{
    std::cerr << "Standard exception: " << e.what() << " ---> Abort!" << std::endl;
    return 4;
}
catch (...)
{
    std::cerr << "Unknown exception thrown! ---> Abort!" << std::endl;
    return 5;
}
