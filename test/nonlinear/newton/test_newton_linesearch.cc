//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>

#include <iostream>
#include <cmath>
#include <cassert>
#include <iomanip>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/istl/bvector.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/nonlinear/newtonsolver_linesearch_cubicbacktracking.hh>

#include "test_newton_common.hh"

int main(int argc, char* argv[])
{
    using namespace Dumux;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);

    // initialize parameters
    Dumux::Parameters::init(argc, argv);

    // use the Newton solver to find a solution to a scalar equation
    using Assembler = MockScalarAssembler;
    using LinearSolver = MockScalarLinearSolver;
    using Solver = NewtonSolver<Assembler, LinearSolver, DefaultPartialReassembler>;

    auto assembler = std::make_shared<Assembler>();
    auto linearSolver = std::make_shared<LinearSolver>();
    auto solver = std::make_shared<Solver>(assembler, linearSolver);

    double initialGuess = 0.1;
    double x = initialGuess;

    std::cout << "Solving: x^2 - 5 = 0" << std::endl;
    solver->solve(x);
    std::cout << "Solution: " << std::setprecision(15) << x
              << ", exact: " << std::sqrt(5.0)
              << ", error: " << std::abs(x-std::sqrt(5.0))/std::sqrt(5.0)*100 << "%" << std::endl;

    if (Dune::FloatCmp::ne(x, std::sqrt(5.0), 1e-13))
        DUNE_THROW(Dune::Exception, "Didn't find correct root: " << std::setprecision(15) << x << ", exact: " << std::sqrt(5.0));

    return 0;

}
