#include <config.h>

#include <iostream>
#include <cmath>
#include <iomanip>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/parallel/mpihelper.hh>
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

class MockScalarVariables
{
public:
    static constexpr int dimension = 1;
    MockScalarVariables() : var_(0.0) {}
    explicit MockScalarVariables(double v) : var_(v) {}
    MockScalarVariables& operator-=(const MockScalarVariables& other) { var_ -= other.var_; return *this; }
    double& operator[] (int i) { return var_; }
    const double& operator[] (int i) const { return var_; }
private:
    double var_;
};

class MockScalarResidual
{
public:
    MockScalarResidual() : res_(0.0) {}
    explicit MockScalarResidual(double r) : res_(r) {}
    MockScalarResidual& operator=(double r) { res_[0] = r; return *this; }
    MockScalarResidual& operator*=(double a) { res_[0] *= a; return *this; }
    MockScalarResidual& operator+=(const MockScalarResidual& other) { res_[0] += other.res_[0]; return *this; }
    MockScalarResidual& operator-=(const MockScalarResidual& other) { res_[0] -= other.res_[0]; return *this; }
    constexpr std::size_t size() const { return 1; }
    MockScalarVariables& operator[] (int i) { return res_; }
    const MockScalarVariables& operator[] (int i) const { return res_; }
    double two_norm2() const { return res_[0]*res_[0]; }
private:
    MockScalarVariables res_;
};

class MockScalarAssembler
{
public:
    using ResidualType = MockScalarResidual;
    using Scalar = double;
    using JacobianMatrix = double;

    void setLinearSystem() {}

    bool isStationaryProblem() { return false; }

    ResidualType prevSol() { return ResidualType(0.0); }

    void resetTimeStep(const ResidualType& sol) {}

    void assembleResidual(const ResidualType& sol)
    {
        res_[0][0] = sol[0][0]*sol[0][0] - 5.0;
    }

    void assembleJacobianAndResidual (const ResidualType& sol)
    {
        assembleResidual(sol);
        jac_ = 2.0*sol[0][0];
    }

    JacobianMatrix& jacobian() { return jac_; }

    ResidualType& residual() { return res_; }

    double residualNorm(const ResidualType& sol)
    {
        assembleResidual(sol);
        return res_[0][0];
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
        x[0][0] = b[0][0]/A;
        return true;
    }
};

} // end namespace Dumux

int main(int argc, char* argv[]) try
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
    MockScalarResidual x(initialGuess);

    std::cout << "Solving: x^2 - 5 = 0" << std::endl;
    solver->solve(x);
    std::cout << "Solution: " << std::setprecision(15) << x[0][0]
              << ", exact: " << std::sqrt(5.0)
              << ", error: " << std::abs(x[0][0]-std::sqrt(5.0))/std::sqrt(5.0)*100 << "%" << std::endl;

    if (Dune::FloatCmp::ne(x[0][0], std::sqrt(5.0), 1e-13))
        DUNE_THROW(Dune::Exception, "Didn't find correct root: " << std::setprecision(15) << x[0][0] << ", exact: " << std::sqrt(5.0));

    return 0;

}
catch (const Dune::Exception& e)
{
    std::cout << e << std::endl;
    return 1;
}
catch (...)
{
    std::cout << "Unknown exception thrown!" << std::endl;
    return 1;
}
