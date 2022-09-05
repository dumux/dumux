#ifndef DUMUX_TEST_TIMESTEPPING_SCALAR_EQUATION_HH
#define DUMUX_TEST_TIMESTEPPING_SCALAR_EQUATION_HH

#include <cmath>
#include <memory>

#include <dune/common/float_cmp.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/istl/multitypeblockvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/io.hh>

#include <dumux/io/format.hh>

#include <dumux/experimental/new_assembly/dumux/common/variables.hh>
#include <dumux/experimental/new_assembly/dumux/common/linearization.hh>

#include <dumux/experimental/new_assembly/dumux/timestepping/multistagevariables.hh>
#include <dumux/experimental/new_assembly/dumux/timestepping/multistagetimestepper.hh>

/*
   This tests the time integration methods by solving the
   linear ODE du/dt = exp(t), where u is the unknown
   and t is the time. We start at t = 0, and thus, the
   exact solution is u_e = exp(t) - 1.
 */

template<typename Dofs>
class TestJacobianOperator : public Dune::LinearOperator<Dofs, Dofs>
{
    using Scalar = Dumux::LinearSystem::ScalarType<Dofs>;

public:
    explicit TestJacobianOperator(Scalar temporalWeight)
    : weight_(temporalWeight)
    {}

    void apply(const Dofs& u, Dofs& r) const override
    {
        r = u;
        r *= weight_;
    }

    void applyscaleadd(Scalar alpha, const Dofs& u, Dofs& r) const override
    {
        auto tmp = r;
        apply(u, tmp);
        tmp *= alpha;
        r += tmp;
    }

    Dune::SolverCategory::Category category() const override
    { return Dune::SolverCategory::Category::sequential; }

private:
    Scalar weight_;
};


// Note: Expects statically sized residual vectors
template<typename V, typename R>
class TestResidualFunction
{
    using Scalar = Dumux::LinearSystem::ScalarType<R>;

public:
    using Domain = Dumux::MultiStageVariables<V>;
    using Range = R;
    using Linearization = Dumux::Linearization<TestJacobianOperator<R>, R>;

    TestResidualFunction()
    : jacOperator_{0.0}
    , residual_{}
    {}

    Range evaluateAt(const Domain& vars)
    {
        Dumux::LinearSystem::fill(residual_, 0.0);
        for (std::size_t stageIdx = 0; stageIdx < vars.numStages(); ++stageIdx)
        {
            addTemporal_(
                residual_,
                vars.stageVariables(stageIdx),
                vars.stageParams().temporalWeight(stageIdx)
            );

            addSpatial_(
                residual_,
                vars.stageVariables(stageIdx),
                vars.stageParams().spatialWeight(stageIdx)
            );
        }
        return residual_;
    }

    Linearization linearizeAt(const Domain& vars)
    {
        evaluateAt(vars);
        const auto curStageIdx = vars.numStages() - 1;
        const auto curTemporalWeight = vars.stageParams().temporalWeight(curStageIdx);
        jacOperator_ = TestJacobianOperator<R>{curTemporalWeight};
        return {jacOperator_, residual_};
    }

private:
    void addTemporal_(Range& residual, const V& v, Scalar weight) const
    { Dumux::LinearSystem::axpy(residual, weight, v.dofs()); }

    void addSpatial_(Range& residual, const V& v, Scalar weight) const
    {
        using std::exp;
        Range tmp;
        Dumux::LinearSystem::fill(tmp, weight*exp(v.timeLevel().current()));
        Dumux::LinearSystem::axpy(residual, -1.0, tmp);
    }

    TestJacobianOperator<R> jacOperator_;
    Range residual_;
};


template<typename Scalar>
Scalar exact(const Scalar t)
{
    using std::exp;
    return exp(t) - 1.0;
}


template<typename Variables>
auto computeError(const Variables& vars)
{
    const auto time = vars.timeLevel().current();
    const auto exactSol = exact(time);

    const auto& dofs = Dumux::Variables::getDofs(vars);
    auto difference = Dumux::Variables::getDofs(vars);
    Dumux::LinearSystem::fill(difference, exactSol);
    Dumux::LinearSystem::subtract(difference, dofs);

    using std::sqrt;
    const auto absErr = sqrt((difference*difference)/dofs.size());
    const auto relErr = absErr/exactSol;
    return std::make_pair(absErr, relErr);
}

template<typename V>
void printVector(const V& v)
{
    Dune::printvector(std::cout, v, "", "");
}

template<typename... Args>
void printVector(const Dune::MultiTypeBlockVector<Args...>& v)
{
    Dune::Hybrid::forEach(std::make_index_sequence<v.size()>{}, [&] (auto i) {
        printVector(v[i]);
    });
}

template<typename Solver, typename Method, typename Scalar>
void testIntegration(std::shared_ptr<Solver> pdeSolver,
                     std::shared_ptr<Method> method,
                     Scalar expectedRelError)
{
    std::cout << "\n-- Integration with " << method->name() << ":\n\n";

    using Variables = typename Solver::Variables;
    using Dofs = Dumux::Variables::DofsType<Variables>;
    using TimeStepper = Dumux::MultiStageTimeStepper<Variables>;

    Dofs dofs;
    dofs = 0.0;
    Variables vars(std::move(dofs));
    TimeStepper timeStepper(pdeSolver, method);

    const Scalar dt = 0.01;
    timeStepper.step(vars, /*time*/0.0, dt);

    const auto [abs, rel] = computeError(vars);
    const auto newTime = vars.timeLevel().current();
    const auto exactSol = exact(newTime);
    if (std::abs(newTime - dt) > 1e-6)
        DUNE_THROW(Dune::InvalidStateException, "Current time level does not match");

    std::cout << "\n"
                << "-- Summary\n"
                << "-- ===========================\n"
                << "-- Discrete solution:\n";
    printVector(Dumux::Variables::getDofs(vars));
    std::cout << "-- Exact solution: " << Dumux::Fmt::format("{:.4e}\n", exactSol)
              << "-- Absolute error: " << Dumux::Fmt::format("{:.4e}\n", abs)
              << "-- Relative error: " << Dumux::Fmt::format("{:.4e}", rel) << std::endl;

    if (Dune::FloatCmp::ne(rel, expectedRelError, 0.01))
        DUNE_THROW(Dune::InvalidStateException,
                    "Detected deviation from reference > 1% for " << method->name());
}

#endif
