// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Convergence test for the complex-valued Helmholtz model against a manufactured solution.
 *        Checks the linear PDE solver and the Newton solver (which has to converge in one step).
 */
#include <config.h>

#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/integrate.hh>

#include <dumux/io/grid/gridmanager_yasp.hh>

#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/pdesolver.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/assembly/fvassembler.hh>

#include "properties.hh"

#ifndef TYPETAG
#define TYPETAG ComplexHelmholtzBox
#endif

#ifndef LINEARSOLVER
#define LINEARSOLVER ILUBiCGSTABIstlSolver
#endif

namespace Dumux {

struct RunResult
{
    double l2Error;
    double h;
};

template<class TypeTag, class GridView>
RunResult solveAndComputeError(const GridView& gridView)
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;

    static_assert(std::is_same_v<typename PrimaryVariables::value_type, std::complex<double>>);
    static_assert(std::is_same_v<GetPropType<TypeTag, Properties::Scalar>, double>);
    static_assert(std::is_same_v<
        typename GetPropType<TypeTag, Properties::JacobianMatrix>::block_type,
        Dune::FieldMatrix<std::complex<double>, 1, 1>
    >, "The default Jacobian block type has to follow the primary variable type");

    auto gridGeometry = std::make_shared<GridGeometry>(gridView);
    auto problem = std::make_shared<Problem>(gridGeometry);
    SolutionVector sol(gridGeometry->numDofs());
    sol = 0.0;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(sol);

    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    using LinearSolver = LINEARSOLVER<LinearSolverTraits<GridGeometry>, LinearAlgebraTraitsFromAssembler<Assembler>>;

    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);
    auto linearSolver = [&]{
        if constexpr (std::is_constructible_v<LinearSolver, const GridView&, const typename GridGeometry::DofMapper&>)
            return std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());
        else
            return std::make_shared<LinearSolver>();
    }();

    LinearPDESolver<Assembler, LinearSolver> linearPDESolver(assembler, linearSolver);
    linearPDESolver.solve(sol);

    // the problem is linear so Newton has to reproduce the solution
    SolutionVector solNewton(sol.size());
    solNewton = 0.0;
    gridVariables->update(solNewton);
    NewtonSolver<Assembler, LinearSolver> newtonSolver(assembler, linearSolver);
    newtonSolver.solve(solNewton);

    // the one-shot linear solve inherits the roundoff error of the numeric Jacobian
    // (of the order of machine precision divided by the numeric epsilon), Newton converges on the residual
    auto diff = sol;
    diff -= solNewton;
    const double relDiff = diff.two_norm()/sol.two_norm();
    std::cout << "Relative difference between Newton and linear solve: " << relDiff << std::endl;
    if (relDiff > 1e-6)
        DUNE_THROW(Dune::Exception, "Newton solution deviates from linear solve: " << relDiff);

    SolutionVector exact(sol.size());
    auto fvGeometry = localView(*gridGeometry);
    for (const auto& element : elements(gridView))
    {
        fvGeometry.bindElement(element);
        for (const auto& scv : scvs(fvGeometry))
            exact[scv.dofIndex()] = PrimaryVariables(problem->exactSolution(scv.dofPosition()));
    }

    const double l2Error = integrateL2Error(*gridGeometry, sol, exact, 4);
    static_assert(std::is_same_v<std::decay_t<decltype(integrateL2Error(*gridGeometry, sol, exact, 4))>, double>,
                  "The L2 error has to be real-valued");

    const double h = std::pow(1.0/gridView.size(0), 1.0/GridView::dimension);
    std::cout << "h = " << h << ", L2 error = " << l2Error << std::endl;
    return { l2Error, h };
}

} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv);

    using TypeTag = Properties::TTag::TYPETAG;
    using Grid = GetPropType<TypeTag, Properties::Grid>;

    GridManager<Grid> gridManager;
    gridManager.init();
    auto& grid = gridManager.grid();

    const int numRefinements = getParam<int>("Grid.Refinements", 2);
    std::vector<double> errors, hs;
    for (int i = 0; i < numRefinements; ++i)
    {
        if (i > 0)
            grid.globalRefine(1);
        const auto result = solveAndComputeError<TypeTag>(grid.leafGridView());
        errors.push_back(result.l2Error);
        hs.push_back(result.h);
    }

    const double maxError = getParam<double>("Problem.MaxL2Error", 5e-3);
    if (errors.back() > maxError)
        DUNE_THROW(Dune::Exception, "L2 error " << errors.back() << " exceeds " << maxError);

    for (std::size_t i = 1; i < errors.size(); ++i)
    {
        const double rate = std::log(errors[i-1]/errors[i])/std::log(hs[i-1]/hs[i]);
        std::cout << "Convergence rate: " << rate << std::endl;
        if (rate < 1.8)
            DUNE_THROW(Dune::Exception, "Convergence rate " << rate << " below 1.8");
    }

    return 0;
}
