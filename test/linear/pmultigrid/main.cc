// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief The Laplace equation with a higher-order Lagrange scheme, solved with the p-multigrid preconditioner
 *        on a sequence of uniformly refined grids, sequentially or in parallel
 *
 * Checks that the L2 error converges with at least the rate Test.MinRate between the two finest grids, that
 * the Krylov method converges to the reduction LinearSolver.ResidualReduction in at most Test.MaxIterations
 * iterations on all grids, and that solving with a stored matrix (setMatrix) gives the same update up to the
 * solver tolerance, also after a solve with another matrix and a reconstruction of the stored solver
 * (setResidualReduction) with the same number of iterations up to one.
 */
#include <config.h>

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/partitionset.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/assembly/assembler.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/grid/gridmanager_alu.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/pmultigridsolvers.hh>
#include <dumux/linear/symmetrize_constraints.hh>

#include "properties.hh"

namespace Dumux {

//! The L2 error of the discrete solution over the interior elements of all processes
template<class GridGeometry, class SolutionVector, class ExactSolution>
double l2Error(const GridGeometry& gridGeometry, const SolutionVector& x, const ExactSolution& exactSolution)
{
    constexpr int dim = GridGeometry::GridView::dimension;
    double error = 0.0;
    auto elemDisc = localView(gridGeometry);
    std::vector<Dune::FieldVector<double, 1>> shapeValues;
    for (const auto& element : elements(gridGeometry.gridView(), Dune::Partitions::interior))
    {
        elemDisc.bind(element);
        const auto geometry = element.geometry();
        const auto& basis = elemDisc.feLocalBasis();
        for (const auto& qp : Dune::QuadratureRules<double, dim>::rule(element.type(), 2*basis.order() + 2))
        {
            basis.evaluateFunction(qp.position(), shapeValues);
            double solution = 0.0;
            for (const auto& localDof : localDofs(elemDisc))
                solution += x[localDof.dofIndex()][0]*shapeValues[localDof.index()][0];
            const double difference = solution - exactSolution(geometry.global(qp.position()));
            error += difference*difference*qp.weight()*geometry.integrationElement(qp.position());
        }
    }
    return std::sqrt(gridGeometry.gridView().comm().sum(error));
}

} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;
    initialize(argc, argv);
    Parameters::init(argc, argv);

    using TypeTag = Properties::TTag::LaplacePMultigrid;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using Assembler = Experimental::Assembler<TypeTag, DiffMethod::numeric>;
    using LinearSolver = SOLVER<LinearSolverTraits<GridGeometry>, LinearAlgebraTraitsFromAssembler<Assembler>>;
    constexpr int dim = Grid::dimension;

    GridManager<Grid> gridManager;
    gridManager.init();
    gridManager.loadBalance();
    const auto& leafGridView = gridManager.grid().leafGridView();
    const bool isRoot = leafGridView.comm().rank() == 0;

    // bisection of simplices halves the element size after dim refinements
    using GridCapabilities = Dune::Capabilities::hasSingleGeometryType<Grid>;
    constexpr bool isSimplexGrid = GridCapabilities::v && GridCapabilities::topologyId == Dune::GeometryTypes::simplex(dim).id();
    const int refinementsPerLevel = isSimplexGrid ? dim : 1;
    const int numLevels = getParam<int>("Test.NumLevels");
    const int maxIterations = getParam<int>("Test.MaxIterations");
    const double minRate = getParam<double>("Test.MinRate");

    bool passed = true;
    double previousError = 0.0;
    for (int level = 0; level < numLevels; ++level)
    {
        if (level > 0)
        {
            gridManager.grid().globalRefine(refinementsPerLevel);
            gridManager.loadBalance();
        }

        auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);
        auto problem = std::make_shared<Problem>(gridGeometry);
        SolutionVector x(gridGeometry->numDofs());
        x = 0.0;
        auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
        gridVariables->init(x);
        auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);
        LinearSolver linearSolver(*gridGeometry);

        // the problem is linear: one Newton step from zero gives the discrete solution;
        // eliminating the constrained columns makes the matrix symmetric for the conjugate gradient method
        assembler->assembleJacobianAndResidual(x);
        auto residual = assembler->residual();
        SolutionVector constrainedRows(x.size());
        constrainedRows = 0.0;
        for (const auto& constraint : problem->constraints())
            constrainedRows[constraint.dofIndex()] = 1.0;
        symmetrizeConstraints(assembler->jacobian(), residual, constrainedRows);

        // the solver prepares matrix and residual in place in a parallel run
        auto matrixAsAssembled = assembler->jacobian();
        const auto residualAsAssembled = residual;
        auto otherMatrix = matrixAsAssembled;
        for (std::size_t i = 0; i < otherMatrix.N(); ++i)
            otherMatrix[i][i] *= 4.0;

        SolutionVector deltaX(x.size());
        deltaX = 0.0;
        const auto result = linearSolver.solve(assembler->jacobian(), deltaX, residual);
        x -= deltaX;

        // the stored-matrix interface gives the same update
        linearSolver.setMatrix(matrixAsAssembled);
        auto residualForReuse = residualAsAssembled;
        SolutionVector deltaXReuse(x.size());
        deltaXReuse = 0.0;
        const auto resultReuse = linearSolver.solve(deltaXReuse, residualForReuse);
        deltaXReuse -= deltaX;
        const double reuseDifference = leafGridView.comm().max(deltaXReuse.infinity_norm())
                                       /leafGridView.comm().max(deltaX.infinity_norm());

        // a solve with another matrix keeps the stored one, also for the reconstructed stored solver
        auto otherResidual = residualAsAssembled;
        SolutionVector otherUpdate(x.size());
        otherUpdate = 0.0;
        linearSolver.solve(otherMatrix, otherUpdate, otherResidual);
        linearSolver.setResidualReduction(getParam<double>("LinearSolver.ResidualReduction"));
        auto residualAfterOther = residualAsAssembled;
        SolutionVector deltaXAfterOther(x.size());
        deltaXAfterOther = 0.0;
        const auto resultAfterOther = linearSolver.solve(deltaXAfterOther, residualAfterOther);
        deltaXAfterOther -= deltaX;
        const double afterOtherDifference = leafGridView.comm().max(deltaXAfterOther.infinity_norm())
                                            /leafGridView.comm().max(deltaX.infinity_norm());
        // sums over more than two processes depend on the order in which their messages arrive
        const int iterationDifference = static_cast<int>(resultAfterOther.iterations)
                                        - static_cast<int>(resultReuse.iterations);
        const bool storedMatrixKept = std::abs(iterationDifference) <= 1 && afterOtherDifference <= 1e-8;

        const double error = l2Error(*gridGeometry, x, Problem::exactSolution);
        const int numElements = leafGridView.comm().sum(static_cast<int>(
            std::distance(elements(leafGridView, Dune::Partitions::interior).begin(),
                          elements(leafGridView, Dune::Partitions::interior).end())));
        const double rate = level > 0 ? std::log2(previousError/error) : 0.0;
        if (isRoot)
            std::cout << std::setw(10) << numElements << " elements: " << std::setw(3) << result.iterations
                      << " iterations (converged " << result.converged << "), L2 error " << std::setprecision(6)
                      << std::scientific << error << ", rate " << std::fixed << std::setprecision(2) << rate << std::endl;

        if (isRoot && reuseDifference > 1e-8)
            std::cout << "stored matrix: " << resultReuse.iterations << " iterations, relative difference of the update "
                      << std::scientific << reuseDifference << std::endl;
        if (isRoot && !storedMatrixKept)
            std::cout << "stored matrix after a solve with another matrix: " << resultAfterOther.iterations
                      << " iterations (before " << resultReuse.iterations << "), relative difference of the update "
                      << std::scientific << afterOtherDifference << std::endl;

        passed = passed && result.converged && result.iterations <= maxIterations && (level < numLevels - 1 || rate >= minRate)
                 && resultReuse.converged && reuseDifference <= 1e-8 && storedMatrixKept;
        previousError = error;
    }

    if (isRoot)
        std::cout << (passed ? "passed" : "failed") << std::endl;
    return passed ? 0 : 1;
}
