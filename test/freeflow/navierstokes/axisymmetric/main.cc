// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Axisymmetric Stokes flow with a manufactured solution: L2 and H1 errors of the velocity
 *        (integrals over the rotated domain)
 */
#include <config.h>

#include <array>
#include <iomanip>
#include <iostream>
#include <memory>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/assembly/assembler.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/io/grid/gridmanager_ug.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <test/freeflow/navierstokes/errors_cvfe.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::TYPETAG;

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv);

    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x;
    problem->applyInitialSolution(x);

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

#if NEW_PROBLEM_INTERFACE
    using Assembler = Experimental::Assembler<TypeTag, DiffMethod::numeric>;
#else
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
#endif
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();
    NewtonSolver<Assembler, LinearSolver> nonLinearSolver(assembler, linearSolver);
    nonLinearSolver.solve(x);

    const auto [volume, errors] = calculateL2AndH1Errors(*problem, *gridVariables, x);
    const auto cells = getParam<std::array<int, 2>>("Grid.Cells");
    std::cout << std::scientific << std::setprecision(10)
              << "[ConvergenceTest] cells " << cells[0]
              << " numDofs " << gridGeometry->numDofs()
              << " volume " << volume
              << " errorL2 " << errors[0]
              << " errorH1 " << errors[1] << std::endl;

    return 0;
}
