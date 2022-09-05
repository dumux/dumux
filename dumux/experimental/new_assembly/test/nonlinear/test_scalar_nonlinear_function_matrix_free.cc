#include <cmath>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/istl/operators.hh>

#include <dumux/common/initialize.hh>
#include <dumux/experimental/new_assembly/dumux/common/linearization.hh>
#include <dumux/experimental/new_assembly/dumux/nonlinear/newton.hh>
#include <dumux/experimental/new_assembly/dumux/linear/dune/solvers.hh>


using Dofs = Dune::FieldVector<double, 1>;

class MyJacobianOperator
: public Dune::LinearOperator<Dofs, Dofs>
{
public:
    explicit MyJacobianOperator(Dofs evalPoint)
    : evalPoint_(evalPoint)
    {}

    void apply(const Dofs& u, Dofs& r) const override
    { r = 2.0*evalPoint_*u; }

    void applyscaleadd(double alpha, const Dofs& u, Dofs& r) const override
    {
        auto tmp = r;
        apply(u, tmp);
        tmp *= alpha;
        r += tmp;
    }

    Dune::SolverCategory::Category category() const override
    { return Dune::SolverCategory::Category::sequential; }

private:
    Dofs evalPoint_;
};

class MyResidualFunction
{
public:
    using Domain = Dofs;
    using Range = Dofs;
    using Linearization = Dumux::Linearization<MyJacobianOperator, Range>;

    const Range& evaluateAt(const Domain& vars)
    {
        residual_ = Dofs{vars*vars - 5.0};
        return residual_;
    }

    Linearization linearizeAt(const Domain& vars)
    {
        evaluateAt(vars);
        operator_ = MyJacobianOperator{vars};
        return {operator_, residual_};
    }

private:
    Range residual_;
    MyJacobianOperator operator_{residual_};
};

int main(int argc, char* argv[])
{
    using namespace Dumux;

    Dumux::initialize(argc, argv);

    NewtonSolver newton{
        std::make_shared<MyResidualFunction>(),
        std::make_shared<DuneCGSolver<Dofs>>()
    };

    Dofs x(0.1);
    if (!newton.solve(x))
        DUNE_THROW(Dune::InvalidStateException, "Newton did not converge");

    using std::abs;
    if (abs(x*x - 5.) > 1e-6)
        DUNE_THROW(Dune::InvalidStateException, "Wrong solution");

    return 0;
}
