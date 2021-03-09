#include <config.h>

#include <iostream>
#include <cmath>
#include <iomanip>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/istl/bvector.hh>
#include <dumux/nonlinear/newtonsolver.hh>

/*

  This test currently solves a scalar non-linear equation using the
  Dumux::NewtonSolver. The Mock classes expose which dependencies the
  current implementation has on different other classes it interacts with.
  Several dependencies seem unnecessary. In particular, the current
  implementation is basically hard-coded to double indexable residual vectors
  with several math operators. A good idea would seem to somehow delegate
  this dependency to something like a linear algebra backend or at least
  the assembler. The assembler requires a lot of interface functions which
  are not always needed. The linear solver interdependency is much better (small).

  This test is to ensure that the dependencies do not grow more in the future.

 */

namespace Dumux {

class MockScalarAssembler
{
public:
    using ResidualType = Dune::BlockVector<double>;
    using Scalar = double;
    using JacobianMatrix = double;

    void setLinearSystem() {}

    bool isStationaryProblem() { return false; }

    ResidualType prevSol() { return ResidualType(0.0); }

    void resetTimeStep(const ResidualType& sol) {}

    void assembleResidual(const ResidualType& sol)
    {
        res_.resize(1);
        res_[0] = sol[0]*sol[0] - 5.0;
    }

    void assembleJacobianAndResidual (const ResidualType& sol)
    {
        assembleResidual(sol);
        jac_ = 2.0*sol[0];
    }

    JacobianMatrix& jacobian() { return jac_; }

    ResidualType& residual() { return res_; }

    double residualNorm(const ResidualType& sol)
    {
        assembleResidual(sol);
        return res_[0];
    }

    void updateGridVariables(const ResidualType& sol) {}

private:
    JacobianMatrix jac_;
    ResidualType res_;
};

class MockScalarLinearSolver
{
public:
    void setResidualReduction(double residualReduction) {}

    template<class Vector>
    bool solve (const double& A, Vector& x, const Vector& b) const
    {
        x[0] = b[0]/A;
        return true;
    }
};

} // end namespace Dumux

int main(int argc, char* argv[])
{
    using namespace Dumux;

    // maybe initialize MPI
    Dune::MPIHelper::instance(argc, argv);

    // initialize  parameters
    // TODO this is necessary because there are some global default used in the Newton solver
    // Do we really need them to be global defaults???
    Parameters::init(argc, argv);

    // use the Newton solver to find a solution to a scalar equation
    using Assembler = MockScalarAssembler;
    using LinearSolver = MockScalarLinearSolver;
    using Solver = NewtonSolver<Assembler, LinearSolver, DefaultPartialReassembler>;

    auto assembler = std::make_shared<Assembler>();
    auto linearSolver = std::make_shared<LinearSolver>();
    auto solver = std::make_shared<Solver>(assembler, linearSolver);

    double initialGuess = 0.1;
    Dune::BlockVector<double> x(1);
    x = initialGuess;

    std::cout << "Solving: x^2 - 5 = 0" << std::endl;
    solver->solve(x);
    std::cout << "Solution: " << std::setprecision(15) << x[0]
              << ", exact: " << std::sqrt(5.0)
              << ", error: " << std::abs(x[0]-std::sqrt(5.0))/std::sqrt(5.0)*100 << "%" << std::endl;

    if (Dune::FloatCmp::ne(x[0], std::sqrt(5.0), 1e-13))
        DUNE_THROW(Dune::Exception, "Didn't find correct root: " << std::setprecision(15) << x[0] << ", exact: " << std::sqrt(5.0));

    return 0;

}
