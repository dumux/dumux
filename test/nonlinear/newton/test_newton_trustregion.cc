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

#include <dumux/nonlinear/newtonsolver_trustregion.hh>

namespace Dumux {

class RosenbrockAssembler
{
public:
    using Scalar = double;
    using ResidualType = Dune::FieldVector<Scalar, 2>;
    using JacobianMatrix = Dune::FieldMatrix<Scalar, 2, 2>;
    using SolutionVector = Dune::FieldVector<Scalar, 2>;
    using Variables = Dune::FieldVector<Scalar, 2>;

    void setLinearSystem() {}

    // Rosenbrock function is f(x) = 100*(x[1] - x[0]^2)^2 + (1 - x[0])^2
    // Residual is the derivative of f'(x) with respect to x[0], x[1], i.e.
    // the function we want to find the root of is f'(x) = 0
    void assembleResidual(const ResidualType& x)
    {
        res_[0] = -400*x[0]*(x[1] - x[0]*x[0]) - 2*(1 - x[0]);
        res_[1] = 200*(x[1] - x[0]*x[0]);
    }

    // the Hessian of the Rosenbrock function
    void assembleJacobianAndResidual (const ResidualType& x)
    {
        assembleResidual(x);
        jac_[0][0] = -400*(x[1] - 3*x[0]*x[0]) + 2;
        jac_[0][1] = -400*x[0];
        jac_[1][0] = -400*x[0];
        jac_[1][1] = 200;
    }

    JacobianMatrix& jacobian() { return jac_; }

    ResidualType& residual() { return res_; }

private:
    JacobianMatrix jac_;
    ResidualType res_;
};

class DirectLinearSolver
{
public:
    void setResidualReduction(double residualReduction) {}

    template <typename Matrix, typename Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b) const
    {
        A.solve(x, b);
        return true;
    }

    template <typename Vector>
    double norm(const Vector& residual) const
    { return residual.two_norm(); }
};

} // end namespace Dumux

int main(int argc, char* argv[])
{
    using namespace Dumux;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);

    // initialize parameters
    Dumux::Parameters::init(argc, argv);

    // use the Newton solver to find a solution to a scalar equation
    using Assembler = RosenbrockAssembler;
    using LinearSolver = DirectLinearSolver;
    using Solver = NewtonSolverTrustRegion<Assembler, LinearSolver>;

    auto assembler = std::make_shared<Assembler>();
    auto linearSolver = std::make_shared<LinearSolver>();
    auto solver = std::make_shared<Solver>(assembler, linearSolver);

    Dune::FieldVector<double, 2> x({1.3, 0.7});
    solver->solve(x);

    std::cout << "\n*************************************\nFound minimum at: (" << x << ")";
    if (!Dune::FloatCmp::eq(x, Dune::FieldVector<double, 2>({1.0, 1.0}), 1e-5))
    {
        const auto error = (Dune::FieldVector<double, 2>({1.0, 1.0}) - x).two_norm();
        std::cout << std::endl;
        DUNE_THROW(Dune::Exception, "Wrong solution: " << x << ", expected (1 1), error: " << error);
        return 1;
    }
    std::cout << " ==> correct\n" << "*************************************" << std::endl;

    return 0;
}
