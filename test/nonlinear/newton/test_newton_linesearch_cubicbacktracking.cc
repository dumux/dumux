//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/nonlinear/newtonsolver_linesearch_cubicbacktracking.hh>

#include "test_newton_common.hh"

namespace Dumux {

double maxError(double x, double exact)
{ return std::abs(x - exact); }

double maxError(Dune::FieldVector<double, 2> x, const Dune::FieldVector<double, 2>& exact)
{ x -= exact; return x.infinity_norm(); }

template<class LinearSolver, class Assembler, class SolutionVector>
void solveAndCheck(const std::string& equation,
                   std::shared_ptr<Assembler> assembler,
                   SolutionVector x,
                   const SolutionVector& exact)
{
    using Solver = NewtonSolverLinesearchCubicbacktracking<Assembler, LinearSolver, DefaultPartialReassembler>;
    Solver solver(assembler, std::make_shared<LinearSolver>());

    std::cout << "Solving: " << equation << std::endl;
    solver.solve(x);

    const auto error = maxError(x, exact);
    std::cout << "Solution: " << std::setprecision(15) << x << ", exact: " << exact << ", error: " << error << std::endl;
    if (!(error < 1e-12))
        DUNE_THROW(Dune::Exception, "Didn't find correct root of " << equation << ": " << x << ", exact: " << exact);
}

} // end namespace Dumux

int main(int argc, char* argv[])
{
    using namespace Dumux;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);

    // initialize parameters
    Dumux::Parameters::init(argc, argv);

    using Vector = Dune::FieldVector<double, 2>;

    // the first Newton step overshoots to x = 25.05
    solveAndCheck<MockScalarLinearSolver>("x^2 - 5 = 0", std::make_shared<MockScalarAssembler>(), 0.1, std::sqrt(5.0));

    // Newton's method without line search alternates between x0 and -x0
    // and the full step does not decrease the residual norm
    const double x0 = 1.391745200270735;
    auto atanAssembler = std::make_shared<MockSystemAssembler>(
        [](double x){ return std::atan(x); },
        [](double x){ return 1.0/(1.0 + x*x); }
    );
    solveAndCheck<MockSystemLinearSolver>("atan(x_i) = 0", atanAssembler, Vector({x0, -x0}), Vector(0.0));

    // the first full Newton step leaves the domain of the logarithm
    auto logAssembler = std::make_shared<MockSystemAssembler>(
        [](double x){ return std::log(x); },
        [](double x){ return 1.0/x; }
    );
    solveAndCheck<MockSystemLinearSolver>("log(x_i) = 0", logAssembler, Vector({3.0, 0.5}), Vector(1.0));

    return 0;
}
