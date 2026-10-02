// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Complex-valued Helmholtz problem in a bulk domain and a facet domain coupled with the box
 *        facet-coupling scheme, solved with the multidomain assembler and a direct solver on the
 *        complex multi-type system. Writes the L2 errors of both domains to a log file.
 */
#include <config.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/facet/gridmanager.hh>
#include <dumux/multidomain/facet/codimonegridadapter.hh>

#include "properties.hh"

namespace Dumux {

//! L2 error of a complex-valued solution relative to the range of the exact solution
template<class GridView, class Problem, class SolutionVector>
double relativeL2Error(const GridView& gridView, const Problem& problem, const SolutionVector& x)
{
    const auto order = getParam<int>("L2Error.QuadratureOrder");
    double error = 0.0;
    double profileMin = std::numeric_limits<double>::max(), profileMax = std::numeric_limits<double>::lowest();
    for (const auto& element : elements(gridView))
    {
        const auto elemSol = elementSolution(element, x, problem.gridGeometry());
        const auto geometry = element.geometry();
        for (auto&& qp : Dune::QuadratureRules<typename GridView::ctype, GridView::dimension>::rule(geometry.type(), order))
        {
            const auto ip = geometry.global(qp.position());
            const auto uh = evalSolution(element, geometry, problem.gridGeometry(), elemSol, ip)[0];
            const auto u = problem.exact(ip);
            profileMin = std::min(profileMin, problem.exactProfile(ip));
            profileMax = std::max(profileMax, problem.exactProfile(ip));
            error += std::norm(uh - u)*qp.weight()*geometry.integrationElement(qp.position());
        }
    }
    return std::sqrt(error)/(std::abs(problem.amplitude())*(profileMax - profileMin));
}

} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv);

    using BulkTypeTag = Properties::TTag::FacetHelmholtzBulk;
    using LowDimTypeTag = Properties::TTag::FacetHelmholtzLowDim;
    using BulkGrid = GetPropType<BulkTypeTag, Properties::Grid>;
    using LowDimGrid = GetPropType<LowDimTypeTag, Properties::Grid>;

    FacetCouplingGridManager<BulkGrid, LowDimGrid> gridManager;
    gridManager.init();
    const auto& bulkGridView = gridManager.template grid<0>().leafGridView();
    const auto& lowDimGridView = gridManager.template grid<1>().leafGridView();

    using BulkGridGeometry = GetPropType<BulkTypeTag, Properties::GridGeometry>;
    using LowDimGridGeometry = GetPropType<LowDimTypeTag, Properties::GridGeometry>;
    auto lowDimGridGeometry = std::make_shared<LowDimGridGeometry>(lowDimGridView);
    // the box facet-coupling grid geometry creates the interior boundary faces from the facet grid
    CodimOneGridAdapter<typename FacetCouplingGridManager<BulkGrid, LowDimGrid>::Embeddings> facetGridAdapter(gridManager.getEmbeddings());
    auto bulkGridGeometry = std::make_shared<BulkGridGeometry>(bulkGridView, lowDimGridView, facetGridAdapter);

    using Traits = Properties::FacetHelmholtzTraits;
    using CouplingManager = typename Traits::CouplingManager;
    auto couplingManager = std::make_shared<CouplingManager>();

    using BulkProblem = GetPropType<BulkTypeTag, Properties::Problem>;
    using LowDimProblem = GetPropType<LowDimTypeTag, Properties::Problem>;
    auto bulkProblem = std::make_shared<BulkProblem>(bulkGridGeometry, couplingManager, "Bulk");
    auto lowDimProblem = std::make_shared<LowDimProblem>(lowDimGridGeometry, couplingManager, "LowDim");

    using MDTraits = typename Traits::MDTraits;
    static_assert(std::is_same_v<typename MDTraits::JacobianMatrix::field_type, std::complex<double>>,
                  "All blocks of the multidomain Jacobian have to be complex-valued");

    using SolutionVector = typename MDTraits::SolutionVector;
    SolutionVector x;
    static const auto bulkId = typename MDTraits::template SubDomain<0>::Index();
    static const auto lowDimId = typename MDTraits::template SubDomain<1>::Index();
    x[bulkId].resize(bulkGridGeometry->numDofs());
    x[lowDimId].resize(lowDimGridGeometry->numDofs());
    x = 0.0;

    auto couplingMapper = std::make_shared<typename Traits::CouplingMapper>();
    couplingMapper->update(*bulkGridGeometry, *lowDimGridGeometry, gridManager.getEmbeddings());
    couplingManager->init(bulkProblem, lowDimProblem, couplingMapper, x);

    using BulkGridVariables = GetPropType<BulkTypeTag, Properties::GridVariables>;
    using LowDimGridVariables = GetPropType<LowDimTypeTag, Properties::GridVariables>;
    auto bulkGridVariables = std::make_shared<BulkGridVariables>(bulkProblem, bulkGridGeometry);
    auto lowDimGridVariables = std::make_shared<LowDimGridVariables>(lowDimProblem, lowDimGridGeometry);
    bulkGridVariables->init(x[bulkId]);
    lowDimGridVariables->init(x[lowDimId]);

    using Assembler = MultiDomainFVAssembler<MDTraits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(bulkProblem, lowDimProblem),
                                                 std::make_tuple(bulkGridGeometry, lowDimGridGeometry),
                                                 std::make_tuple(bulkGridVariables, lowDimGridVariables),
                                                 couplingManager);

    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();
    MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager> newtonSolver(assembler, linearSolver, couplingManager);
    newtonSolver.solve(x);

    const auto bulkError = relativeL2Error(bulkGridView, *bulkProblem, x[bulkId]);
    const auto lowDimError = relativeL2Error(lowDimGridView, *lowDimProblem, x[lowDimId]);
    const auto epsilon = getParam<double>("LowDim.SpatialParams.Aperture")*getParam<int>("Grid.NumElemsPerSide");
    std::cout << "Matrix - epsilon/error: " << epsilon << ", " << bulkError << std::endl;
    std::cout << "Fracture - epsilon/error: " << epsilon << ", " << lowDimError << std::endl;

    std::ofstream file(getParam<std::string>("Problem.OutputFileName"), std::ios::app);
    file << epsilon << "," << bulkError << "," << epsilon << "," << lowDimError << "\n";
    return 0;
}
