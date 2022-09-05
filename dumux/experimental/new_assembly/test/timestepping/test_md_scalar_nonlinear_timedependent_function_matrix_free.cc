#include <memory>

#include <dune/common/fvector.hh>
#include <dune/istl/multitypeblockvector.hh>

#include <dumux/common/initialize.hh>

#include <dumux/experimental/new_assembly/dumux/common/variables.hh>
#include <dumux/experimental/new_assembly/dumux/nonlinear/newton.hh>
#include <dumux/experimental/new_assembly/dumux/linear/dune/solvers.hh>

#include <dumux/experimental/new_assembly/dumux/timestepping/multistagemethods.hh>

#include "scalar_equation.hh"

using MyScalar = double;
using MyDofs = Dune::MultiTypeBlockVector<Dune::FieldVector<MyScalar, 1>,
                                          Dune::FieldVector<MyScalar, 1>>;
using MyVariables = Dumux::DefaultVariables<MyDofs>;

int main(int argc, char* argv[])
{
    Dumux::initialize(argc, argv);

    using Function = TestResidualFunction<MyVariables, MyDofs>;
    using LinearSolver = Dumux::DuneCGSolver<MyDofs>;
    using NewtonSolver = Dumux::NewtonSolver<Function, LinearSolver>;
    auto newton = std::make_shared<NewtonSolver>(
        std::make_shared<Function>(),
        std::make_shared<LinearSolver>()
    );

    using namespace Dumux::MultiStage;
    testIntegration(newton, std::make_shared<ExplicitEuler<MyScalar>>(), 4.9917e-03);
    testIntegration(newton, std::make_shared<ImplicitEuler<MyScalar>>(), 5.0083e-03);
    testIntegration(newton, std::make_shared<Theta<MyScalar>>(0.5), 8.3333e-06);
    testIntegration(newton, std::make_shared<RungeKuttaExplicitFourthOrder<MyScalar>>(), 3.4829e-12);

    return 0;
}
