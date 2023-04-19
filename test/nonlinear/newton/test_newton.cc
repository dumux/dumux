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
#include <dumux/nonlinear/newtonsolver.hh>

/*
  This test currently solves a scalar non-linear equation using the
  Dumux::NewtonSolver. The Mock classes expose which dependencies the
  current implementation has on different other classes it interacts with.
  This test is to ensure that the dependencies do not grow more in the future.
  The Mock classes can be step-wise reduced in complexity once the Newton
  implementation required a smaller interface from assembler and solver.
 */

namespace Dumux {

class MockScalarAssembler
{
public:
    using Scalar = double;
    using ResidualType = Scalar;
    using JacobianMatrix = Scalar;
    using SolutionVector = Scalar;
    using Variables = Scalar;

    void setLinearSystem() {}

    void assembleResidual(const ResidualType& sol)
    {
        res_ = sol*sol - 5.0;
    }

    void assembleJacobianAndResidual (const ResidualType& sol)
    {
        assembleResidual(sol);
        jac_ = 2.0*sol;
    }

    JacobianMatrix& jacobian() { return jac_; }

    ResidualType& residual() { return res_; }

private:
    JacobianMatrix jac_;
    ResidualType res_;
};

class MockScalarLinearSolver
{
public:
    void setResidualReduction(double residualReduction) {}

    bool solve(const double& A, double& x, const double& b) const
    {
        x = b/A;
        return true;
    }

    double norm(const double& residual) const
    {
        using std::abs;
        return abs(residual);
    }
};

} // end namespace Dumux

int main(int argc, char* argv[])
{
    using namespace Dumux;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);

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
