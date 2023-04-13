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

#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/linear/seqsolverbackend.hh>

#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>

#include <dumux/multidomain/facet/gridmanager.hh>
#include <dumux/multidomain/facet/codimonegridadapter.hh>

#include <dumux/io/vtkoutputmodule.hh>

#include "properties.hh"

namespace Dumux {

//! Container to store the L2 errors
template< class Scalar >
struct L2NormData
{
    Scalar epsilon; // The ratio of aperture/dy
    Scalar norm;    // The L2 norm of the error
};

//! Computes the L2 norm of the error of a sub-problem
template< class GridView, class Problem, class SolutionVector >
L2NormData< typename SolutionVector::field_type >
computeL2Norm(const GridView& gridView,
              const Problem& problem,
              const SolutionVector& x)
{
    // container to store results in
    using Scalar = typename SolutionVector::field_type;
    L2NormData< Scalar > l2Data;

    // keep track of min/max solution values
    Scalar uMin = std::numeric_limits<Scalar>::max();
    Scalar uMax = std::numeric_limits<Scalar>::min();

    // some necessary parameters
    const auto a = Dumux::getParam<Scalar>("LowDim.SpatialParams.Aperture");
    const auto numE = Dumux::getParam<int>("Grid.NumElemsPerSide");
    const auto order = Dumux::getParam<int>("L2Error.QuadratureOrder");

    // Compute epsilon (a/dx = a/(1/numE) = a*numE)
    l2Data.epsilon = a*numE;
    l2Data.norm = 0.0;

    // integrate error in each cell
    for (const auto& element : elements(gridView))
    {
        // make discrete element solution
        const auto elemSol = elementSolution(element, x, problem.gridGeometry());

        // integrate the pressure error over the element
        const auto eg = element.geometry();
        const auto& rule = Dune::QuadratureRules<Scalar, GridView::dimension>::rule(eg.type(), order);
        for (auto&& qp : rule)
        {
            // integration point in global coordinates
            const auto ip = eg.global(qp.position());

            // exact and discrete solution at integration point
            const auto u = problem.exact(ip);
            const auto uh = evalSolution(element, eg, elemSol, ip);

            using std::min;
            using std::max;
            uMin = min(u, uMin);
            uMax = max(u, uMax);
            l2Data.norm += (uh - u)*(uh - u)*qp.weight()*eg.integrationElement(qp.position());
        }
    }

    // take the root of the norm and scale it
    using std::sqrt;
    l2Data.norm = sqrt(l2Data.norm)/(uMax-uMin);

    // return result
    return l2Data;
}

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
    using BulkProblemTypeTag = Properties::TTag::BULKTYPETAG;
    using LowDimProblemTypeTag = Properties::TTag::LOWDIMTYPETAG;
    using BulkGrid = GetPropType<BulkProblemTypeTag, Properties::Grid>;
    using LowDimGrid = GetPropType<LowDimProblemTypeTag, Properties::Grid>;

    using GridManager = FacetCouplingGridManager<BulkGrid, LowDimGrid>;
    GridManager gridManager;
    gridManager.init();
    gridManager.loadBalance();

    ////////////////////////////////////////////////////////////
    // run stationary, non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid views
    const auto& bulkGridView = gridManager.template grid<0>().leafGridView();
    const auto& lowDimGridView = gridManager.template grid<1>().leafGridView();

    // create the finite volume grid geometries
    using BulkFVGridGeometry = GetPropType<BulkProblemTypeTag, Properties::GridGeometry>;
    using LowDimFVGridGeometry = GetPropType<LowDimProblemTypeTag, Properties::GridGeometry>;
    auto lowDimFvGridGeometry = std::make_shared<LowDimFVGridGeometry>(lowDimGridView);
    auto bulkFvGridGeometry  = makeBulkFVGridGeometry<BulkFVGridGeometry>(gridManager, bulkGridView, lowDimGridView);

    // the coupling manager
    using TestTraits = Properties::TestTraits<BulkProblemTypeTag, LowDimProblemTypeTag>;
    using CouplingManager = typename TestTraits::CouplingManager;
    auto couplingManager = std::make_shared<CouplingManager>();

    // the problems (boundary conditions)
    using BulkProblem = GetPropType<BulkProblemTypeTag, Properties::Problem>;
    using LowDimProblem = GetPropType<LowDimProblemTypeTag, Properties::Problem>;
    auto bulkSpatialParams = std::make_shared<typename BulkProblem::SpatialParams>(bulkFvGridGeometry, "Bulk");
    auto bulkProblem = std::make_shared<BulkProblem>(bulkFvGridGeometry, bulkSpatialParams, couplingManager, "Bulk");
    auto lowDimSpatialParams = std::make_shared<typename LowDimProblem::SpatialParams>(lowDimFvGridGeometry, "LowDim");
    auto lowDimProblem = std::make_shared<LowDimProblem>(lowDimFvGridGeometry, lowDimSpatialParams, couplingManager, "LowDim");

    // the solution vector
    using MDTraits = typename TestTraits::MDTraits;
    using SolutionVector = typename MDTraits::SolutionVector;
    SolutionVector x;

    static const auto bulkId = typename MDTraits::template SubDomain<0>::Index();
    static const auto lowDimId = typename MDTraits::template SubDomain<1>::Index();
    x[bulkId].resize(bulkFvGridGeometry->numDofs());
    x[lowDimId].resize(lowDimFvGridGeometry->numDofs());
    bulkProblem->applyInitialSolution(x[bulkId]);
    lowDimProblem->applyInitialSolution(x[lowDimId]);

    // the coupling mapper
    auto couplingMapper = std::make_shared<typename TestTraits::CouplingMapper>();
    couplingMapper->update(*bulkFvGridGeometry, *lowDimFvGridGeometry, gridManager.getEmbeddings());

    // initialize the coupling manager
    couplingManager->init(bulkProblem, lowDimProblem, couplingMapper, x);

    // the grid variables
    using BulkGridVariables = GetPropType<BulkProblemTypeTag, Properties::GridVariables>;
    using LowDimGridVariables = GetPropType<LowDimProblemTypeTag, Properties::GridVariables>;
    auto bulkGridVariables = std::make_shared<BulkGridVariables>(bulkProblem, bulkFvGridGeometry);
    auto lowDimGridVariables = std::make_shared<LowDimGridVariables>(lowDimProblem, lowDimFvGridGeometry);
    bulkGridVariables->init(x[bulkId]);
    lowDimGridVariables->init(x[lowDimId]);

    // initialize the vtk output modulell
    const auto bulkDM = BulkFVGridGeometry::discMethod == DiscretizationMethods::box ? Dune::VTK::nonconforming : Dune::VTK::conforming;
    using BulkSolutionVector = std::decay_t<decltype(x[bulkId])>;
    using LowDimSolutionVector = std::decay_t<decltype(x[lowDimId])>;
    VtkOutputModule<BulkGridVariables, BulkSolutionVector> bulkVtkWriter(*bulkGridVariables, x[bulkId], bulkProblem->name(), "Bulk", bulkDM);
    VtkOutputModule<LowDimGridVariables, LowDimSolutionVector> lowDimVtkWriter(*lowDimGridVariables, x[lowDimId], lowDimProblem->name(), "LowDim");

    // container for the output of the exact solutions
    std::vector<GetPropType<BulkProblemTypeTag, Properties::Scalar>> bulkExact;
    std::vector<GetPropType<LowDimProblemTypeTag, Properties::Scalar>> lowDimExact;

    // Add model specific output fields
    const bool writeVTK = getParam<bool>("Output.EnableVTK");
    if (writeVTK)
    {
        using BulkIOFields = GetPropType<BulkProblemTypeTag, Properties::IOFields>;
        using LowDimIOFields = GetPropType<LowDimProblemTypeTag, Properties::IOFields>;
        BulkIOFields::initOutputModule(bulkVtkWriter);
        LowDimIOFields::initOutputModule(lowDimVtkWriter);

        bulkExact.resize(bulkFvGridGeometry->numDofs());
        lowDimExact.resize(lowDimFvGridGeometry->numDofs());

        for (const auto& element : elements(bulkFvGridGeometry->gridView()))
        {
            if (BulkFVGridGeometry::discMethod == DiscretizationMethods::box)
                for (int i = 0; i < element.geometry().corners(); ++i)
                    bulkExact[ bulkFvGridGeometry->vertexMapper().subIndex(element, i, BulkGrid::dimension) ]
                            = bulkProblem->exact( element.template subEntity<BulkGrid::dimension>(i).geometry().center() );
            else
                bulkExact[ bulkFvGridGeometry->elementMapper().index(element) ] = bulkProblem->exact( element.geometry().center() );
        }

        for (const auto& element : elements(lowDimFvGridGeometry->gridView()))
        {
            if (LowDimFVGridGeometry::discMethod == DiscretizationMethods::box)
                for (int i = 0; i < element.geometry().corners(); ++i)
                    lowDimExact[ lowDimFvGridGeometry->vertexMapper().subIndex(element, i, LowDimGrid::dimension) ]
                            = lowDimProblem->exact( element.template subEntity<LowDimGrid::dimension>(i).geometry().center() );
            else
                lowDimExact[ lowDimFvGridGeometry->elementMapper().index(element) ] = lowDimProblem->exact( element.geometry().center() );
        }

        bulkVtkWriter.addField(bulkExact, "pressure_exact");
        lowDimVtkWriter.addField(lowDimExact, "pressure_exact");

        // write initial solution
        bulkVtkWriter.write(0.0);
        lowDimVtkWriter.write(0.0);
    }

    // the assembler
    using Assembler = MultiDomainFVAssembler<MDTraits, CouplingManager, DiffMethod::numeric, /*implicit?*/true>;
    auto assembler = std::make_shared<Assembler>( std::make_tuple(bulkProblem, lowDimProblem),
                                                  std::make_tuple(bulkFvGridGeometry, lowDimFvGridGeometry),
                                                  std::make_tuple(bulkGridVariables, lowDimGridVariables),
                                                  couplingManager);

    // the linear solver
    using LinearSolver = BlockDiagAMGBiCGSTABSolver;
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
    if (writeVTK)
    {
        bulkVtkWriter.write(1.0);
        lowDimVtkWriter.write(1.0);
    }

    // norm of the errors in the two domains
    const auto m = computeL2Norm(bulkFvGridGeometry->gridView(), *bulkProblem, x[bulkId]);
    const auto f = computeL2Norm(lowDimFvGridGeometry->gridView(), *lowDimProblem, x[lowDimId]);

    // output to terminal and file
    std::cout << "Matrix - epsilon/error: " << m.epsilon << ", " << m.norm << std::endl;
    std::cout << "Fracture - epsilon/error: " << f.epsilon << ", " << f.norm << std::endl;

    const auto outputFileName = getParam<std::string>("Problem.OutputFileName");
    std::ofstream file;
    file.open(outputFileName, std::ios::app);
    file << m.epsilon << "," << m.norm << "," << f.epsilon << "," << f.norm << "\n";
    file.close();

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/false);

    return 0;
}
