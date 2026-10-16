// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Momentum balance of the Donea benchmark (Stokes flow, prescribed pressure) with a higher-order scheme:
 *        BiCGSTAB preconditioned by AMG and by the p-multigrid preconditioner on a sequence of refined grids,
 *        sequentially or in parallel
 *
 * Both solvers solve the same assembled system. Checks that p-multigrid converges to the reduction
 * LinearSolver.ResidualReduction in at most Test.MaxIterations iterations on all grids, that its solution agrees
 * with that of AMG to Test.SolutionTolerance relative to the solution where AMG converges, and that the L2 error
 * of the velocity converges with at least the rate Test.MinRate between the two finest grids.
 */
#include <config.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/istl/istlexception.hh>
#include <dune/istl/solver.hh>
#include <dune/common/timer.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/common/partitionset.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/assembly/assembler.hh>
#include <dumux/io/grid/gridmanager_alu.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/parallelhelpers.hh>
#include <dumux/linear/pmultigridsolvers.hh>

#include "properties_momentum.hh"

namespace Dumux {

//! The L2 error of the velocity over the interior elements of all processes
template<class GridGeometry, class Problem, class SolutionVector>
double velocityL2Error(const GridGeometry& gridGeometry, const Problem& problem, const SolutionVector& x)
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
            Dune::FieldVector<double, dim> difference(0.0);
            for (const auto& localDof : localDofs(elemDisc))
                difference.axpy(shapeValues[localDof.index()][0], x[localDof.dofIndex()]);
            difference -= problem.analyticalSolution(geometry.global(qp.position()));
            error += difference.two_norm2()*qp.weight()*geometry.integrationElement(qp.position());
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

    using TypeTag = Properties::TTag::TYPETAG_MOMENTUM;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using Assembler = Experimental::Assembler<TypeTag, DiffMethod::numeric>;
    using LSTraits = LinearSolverTraits<GridGeometry>;
    using LATraits = LinearAlgebraTraitsFromAssembler<Assembler>;
    using AMGSolver = AMGBiCGSTABIstlSolver<LSTraits, LATraits>;
    using PMultigridSolver = PMultigridBiCGSTABIstlSolver<LSTraits, LATraits>;

    GridManager<Grid> gridManager;
    gridManager.init();
    gridManager.loadBalance();
    const auto& leafGridView = gridManager.grid().leafGridView();
    const bool isRoot = leafGridView.comm().rank() == 0;

    const int numLevels = getParam<int>("Test.NumLevels");
    const int maxIterations = getParam<int>("Test.MaxIterations");
    const double solutionTolerance = getParam<double>("Test.SolutionTolerance");
    const double minRate = getParam<double>("Test.MinRate");

    bool passed = true;
    double previousError = 0.0;
    for (int level = 0; level < numLevels; ++level)
    {
        if (level > 0)
        {
            gridManager.grid().globalRefine(1);
            gridManager.loadBalance();
        }

        auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);
        auto problem = std::make_shared<Problem>(gridGeometry);
        SolutionVector x(gridGeometry->numDofs());
        x = 0.0;
        auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
        gridVariables->init(x);
        auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

        // the problem is linear: one Newton step from zero gives the discrete solution;
        // the solvers prepare matrix and residual in place in a parallel run, so each gets a copy
        assembler->assembleJacobianAndResidual(x);
        const auto solve = [&](auto& linearSolver, SolutionVector& update)
        {
            auto A = assembler->jacobian();
            auto b = assembler->residual();
            update.resize(x.size());
            update = 0.0;
            Dune::Timer timer;
            Dune::InverseOperatorResult result;
            // a breakdown of the Krylov method counts as not converged
            try { result = linearSolver.solve(A, update, b); }
            catch (const Dune::ISTLError&) { result.converged = false; }
            return std::make_pair(result, leafGridView.comm().max(timer.elapsed()));
        };

        AMGSolver amgSolver(leafGridView, LSTraits::dofMapper(*gridGeometry));
        SolutionVector amgUpdate;
        const auto [amgResult, amgSeconds] = solve(amgSolver, amgUpdate);

        PMultigridSolver pmgSolver(*gridGeometry);
        SolutionVector pmgUpdate;
        const auto [pmgResult, pmgSeconds] = solve(pmgSolver, pmgUpdate);

        // the values at dofs of ghost entities are not part of the solution
        const ParallelISTLHelper<LSTraits> parallelHelper(leafGridView, LSTraits::dofMapper(*gridGeometry));
        double maxDifference = 0.0, maxUpdate = 0.0;
        for (std::size_t i = 0; i < x.size(); ++i)
        {
            if (parallelHelper.isGhost(i))
                continue;
            auto difference = pmgUpdate[i];
            difference -= amgUpdate[i];
            maxDifference = std::max(maxDifference, difference.infinity_norm());
            maxUpdate = std::max(maxUpdate, amgUpdate[i].infinity_norm());
        }
        const double relativeDifference = leafGridView.comm().max(maxDifference)/leafGridView.comm().max(maxUpdate);

        auto amgSolution = x; amgSolution -= amgUpdate;
        auto pmgSolution = x; pmgSolution -= pmgUpdate;
        const double amgError = velocityL2Error(*gridGeometry, *problem, amgSolution);
        const double pmgError = velocityL2Error(*gridGeometry, *problem, pmgSolution);
        const double rate = level > 0 ? std::log2(previousError/pmgError) : 0.0;

        const auto numDofs = leafGridView.comm().sum(gridGeometry->numDofs());
        if (isRoot)
            std::cout << std::setw(9) << numDofs << " dofs (sum over processes): AMG " << std::setw(3) << amgResult.iterations
                      << " iterations (converged " << amgResult.converged << "), " << std::setprecision(3) << amgSeconds
                      << " s, velocity error " << std::scientific << amgError << std::defaultfloat << "; p-multigrid "
                      << std::setw(3) << pmgResult.iterations << " iterations (converged " << pmgResult.converged << "), "
                      << pmgSeconds << " s, velocity error " << std::scientific << pmgError << ", rate " << std::fixed
                      << std::setprecision(2) << rate << "; relative difference " << std::scientific << relativeDifference
                      << std::defaultfloat << std::endl;

        passed = passed && pmgResult.converged && pmgResult.iterations <= maxIterations
                 && (!amgResult.converged || relativeDifference <= solutionTolerance)
                 && (level < numLevels - 1 || rate >= minRate);
        previousError = pmgError;
    }

    if (isRoot)
        std::cout << (passed ? "passed" : "failed") << std::endl;
    return passed ? 0 : 1;
}
