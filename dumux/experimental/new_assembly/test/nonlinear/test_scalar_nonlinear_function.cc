#include <config.h>

#include <cmath>
#include <memory>

#if HAVE_EIGEN3
#include <Eigen/Core>
#endif

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/exceptions.hh>
#include <dune/istl/operators.hh>

#include <dumux/common/initialize.hh>
#include <dumux/experimental/new_assembly/dumux/common/linearization.hh>
#include <dumux/experimental/new_assembly/dumux/nonlinear/newton.hh>
#include <dumux/experimental/new_assembly/dumux/linear/system.hh>

#include <dumux/experimental/new_assembly/dumux/linear/dune/solvers.hh>
#include <dumux/experimental/new_assembly/dumux/linear/eigen/solvers.hh>

template<typename Vector,
         typename Matrix,
         std::constructible_from<Matrix&> Operator>
class MyResidualFunction
{
public:
    using Domain = Vector;
    using Range = Vector;
    using Linearization = Dumux::Linearization<Operator, Vector>;

    const Range& evaluateAt(const Domain& vars)
    {
        Dumux::LinearSystem::fill(residual_, vars*vars - 5.0);
        return residual_;
    }

    Linearization linearizeAt(const Domain& vars)
    {
        evaluateAt(vars);
        matrix_ = Matrix{2.0*vars};
        jacOperator_ = std::make_unique<Operator>(matrix_);
        return {*jacOperator_, residual_};
    }

private:
    Range residual_;
    Matrix matrix_;
    std::unique_ptr<Operator> jacOperator_{nullptr};
};

template<typename Dofs>
bool isCorrectResult(const Dofs& x)
{
    using std::abs;
    return abs(x*x - 5.0) < 1e-6;
}

template<typename Dofs, typename... Args>
void solveAndCheck(Dofs& x, Dumux::NewtonSolver<Args...>& newtonSolver)
{
    if (!newtonSolver.solve(x))
        DUNE_THROW(Dune::InvalidStateException, "Newton did not converge");

    if (!isCorrectResult(x))
        DUNE_THROW(Dune::InvalidStateException, "Incorrect result");
}

inline void solveWithDune()
{
    std::cout << "Solving with Dune" << std::endl;

    using Dofs = Dune::FieldVector<double, 1>;
    using Matrix = Dune::FieldMatrix<double, 1, 1>;
    using Operator = Dune::MatrixAdapter<Matrix, Dofs, Dofs>;
    using LinearSolver = Dumux::DuneCGSolver<Dofs>;
    using Function = MyResidualFunction<Dofs, Matrix, Operator>;

    Dumux::NewtonSolver newton{
        std::make_shared<Function>(),
        std::make_shared<LinearSolver>()
    };

    Dofs x(0.1);
    solveAndCheck(x, newton);
}

inline void solveWithEigen()
{
#if HAVE_EIGEN3
    std::cout << "Solving with Eigen" << std::endl;

    using Dofs = Eigen::Matrix<double, 1, 1>;
    using Matrix = Eigen::Matrix<double, 1, 1>;
    using LinearSolver = Dumux::EigenCGSolver<Dofs, Matrix>;
    using Function = MyResidualFunction<Dofs, Matrix, typename LinearSolver::LinearOperator>;

    Dumux::NewtonSolver newton{
        std::make_shared<Function>(),
        std::make_shared<LinearSolver>()
    };

    Dofs x(0.1);
    solveAndCheck(x, newton);
#else
    std::cout << "Eigen not found, skipping..." << std::endl;
#endif
}

int main(int argc, char* argv[])
{
    Dumux::initialize(argc, argv);

    solveWithDune();
    solveWithEigen();

    return 0;
}
