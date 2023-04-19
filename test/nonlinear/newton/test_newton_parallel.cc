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

#include <dumux/common/timeloop.hh>
#include <dumux/common/initialize.hh>
#include <dumux/nonlinear/newtonsolver.hh>

/*
  Parallel Newton test. If the solver fails, we reduce the time step size.
  In parallel this can lead to dead locks if we don't correctly handle exceptions.
  We can't avoid all such situation if code in Dune throw and communicates after.
  However, we avoid it in situation where we can recover.
 */

namespace Dumux {

class MockScalarAssembler
{
public:
    using Scalar = double;
    using ResidualType = Scalar;
    using JacobianMatrix = Scalar;
    using SolutionVector = Scalar;

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

    // the following four methods are obsolete for the new assembly
    // chosen if the assembler exports 'Variables' (see test_newton.cc)
    // {
    void updateGridVariables(const ResidualType& sol) {}

    void resetTimeStep(const ResidualType& sol) {}

    ResidualType& prevSol() { return res_; }

    // we fake this for this test
    bool isStationaryProblem() { return false; }
    // }

private:
    JacobianMatrix jac_;
    ResidualType res_;
};

class MockScalarLinearSolver
{
public:
    MockScalarLinearSolver(int rank, std::shared_ptr<TimeLoop<double>> timeLoop)
    : rank_(rank), timeLoop_(timeLoop) {}

    void setResidualReduction(double residualReduction) {}

    template<class Vector>
    bool solve(const double& A, Vector& x, const Vector& b) const
    {
        if (rank_ == 1 && timeLoop_->timeStepSize() > 0.4)
            DUNE_THROW(Dune::Exception, "Assembly failed (for testing) on process " << rank_);

        x = b/A;
        return true;
    }

    double norm(const double& residual) const
    {
        using std::abs;
        return abs(residual);
    }
private:
    int rank_;
    std::shared_ptr<TimeLoop<double>> timeLoop_;
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

    const int rank = Dune::MPIHelper::instance().rank();
    auto timeLoop = std::make_shared<TimeLoop<double>>(0.0, 1.0, 1.0);

    auto assembler = std::make_shared<Assembler>();
    auto linearSolver = std::make_shared<LinearSolver>(rank, timeLoop);
    auto solver = std::make_shared<Solver>(assembler, linearSolver);

    double initialGuess = 0.1;
    double x = initialGuess;

    std::cout << "Solving: x^2 - 5 = 0" << std::endl;
    solver->solve(x, *timeLoop);

    if (rank == 0)
    {
        // check if there have been two time step reductions
        if (Dune::FloatCmp::ne(timeLoop->timeStepSize(), 0.25, 1e-8))
            DUNE_THROW(Dune::Exception, "Time step reduction and recovery did not work");

        std::cout << "Solution: " << std::setprecision(15) << x
                << ", exact: " << std::sqrt(5.0)
                << ", error: " << std::abs(x-std::sqrt(5.0))/std::sqrt(5.0)*100 << "%" << std::endl;

        if (Dune::FloatCmp::ne(x, std::sqrt(5.0), 1e-13))
            DUNE_THROW(Dune::Exception, "Didn't find correct root: " << std::setprecision(15) << x << ", exact: " << std::sqrt(5.0));
    }

    return 0;
}
