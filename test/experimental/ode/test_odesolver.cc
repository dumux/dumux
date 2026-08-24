// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <string>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/experimental/ode/odesolver.hh>

namespace Dumux {

// ODE of form: \f$ \frac{\mathrm{d}u(\xi)}{\mathrm{d}\xi} = \mathrm{e}^\xi \f$
class ExponentialRhsODE
{
public:
    using Scalar = double;
    using SolutionVector = Scalar;
    using ResidualType = Scalar;
    using JacobianMatrix = Scalar;
    using Variables = Experimental::ODEVariables<SolutionVector>;

    void rhs(const Variables& vars, ResidualType& rhs) const
    {
        using std::exp;
        rhs = exp(vars.independentVariableLevel().current());
    }
};

// ODE of form: \f$ \frac{\mathrm{d}u(\xi)}{\mathrm{d}\xi} = \lambda*u(\xi) + \mathrm{e}^\xi \f$
class LinearScalarODE
{
public:
    using Scalar = double;
    using SolutionVector = Scalar;
    using ResidualType = Scalar;
    using JacobianMatrix = Scalar;
    using Variables = Experimental::ODEVariables<SolutionVector>;

    explicit LinearScalarODE(const Scalar lambda)
    : lambda_(lambda)
    {}

    void rhs(const Variables& vars, ResidualType& rhs) const
    {
        using std::exp;
        rhs = lambda_*vars.dofs() + exp(vars.independentVariableLevel().current());
    }

    void rhsJacobian(const Variables&, JacobianMatrix& jacobian) const
    { jacobian = lambda_; }

private:
    Scalar lambda_;
};

/* ODE of form:
 * \f[
 * \begin{aligned}
 * \frac{\mathrm{d}u_1(\xi)}{\mathrm{d}\xi} &= u_2(\xi), \\
 * \frac{\mathrm{d}u_2(\xi)}{\mathrm{d}\xi} &= -u_1(\xi)
 * \end{aligned}
 * \f]
 */
class HarmonicOscillatorODE
{
public:
    using Scalar = double;
    using SolutionVector = Dune::FieldVector<Scalar, 2>;
    using ResidualType = SolutionVector;
    using JacobianMatrix = Dune::FieldMatrix<Scalar, 2, 2>;
    using Variables = Experimental::ODEVariables<SolutionVector>;

    void rhs(const Variables& vars, ResidualType& rhs) const
    {
        rhs[0] = vars.dofs()[1];
        rhs[1] = -vars.dofs()[0];
    }

    void rhsJacobian(const Variables&, JacobianMatrix& jacobian) const
    {
        jacobian = 0.0;
        jacobian[0][1] = 1.0;
        jacobian[1][0] = -1.0;
    }
};

// ODE of form: \f$ \frac{\mathrm{d} (2u(\xi))}{\mathrm{d}\xi} = 1.0 \f$
class ScaledDerivativeQuantityODE
{
public:
    using Scalar = double;
    using SolutionVector = Scalar;
    using ResidualType = Scalar;
    using JacobianMatrix = Scalar;
    using Variables = Experimental::ODEVariables<SolutionVector>;

    void derivativeQuantity(const Variables& vars, ResidualType& value) const
    { value = 2.0*vars.dofs(); }

    void derivativeQuantityJacobian(const Variables&, JacobianMatrix& jacobian) const
    { jacobian = 2.0; }

    void rhs(const Variables&, ResidualType& rhs) const
    { rhs = 1.0; }

    void rhsJacobian(const Variables&, JacobianMatrix& jacobian) const
    { jacobian = 0.0; }
};

} // end namespace Dumux

int main(int argc, char* argv[])
{
    using namespace Dumux;
    using Scalar = double;

    Dumux::initialize(argc, argv);
    Dumux::Parameters::init(argc, argv);

    const auto expectNear = [] (const auto value,
                                const auto reference,
                                const auto tolerance,
                                const std::string& message)
    {
        using std::abs;
        if (abs(value - reference) > tolerance)
            DUNE_THROW(Dune::InvalidStateException,
                       message << " (value: " << value << ", reference: " << reference << ")");
    };

    /*
     * Solves the ODE: \f$ \frac{\mathrm{d}u(\xi)}{\mathrm{d}\xi} = \mathrm{e}^\xi \f$
     * on \f$ \xi \in [0,1] \f$ with \f$ u(0)=0 \f$
     * and a step size of \f$ \Delta \xi = 0.01 \f$
     *
     * note: `ExponentialRhsODE` does not need to define `rhsJacobian()` since an explicit solving strategy is used
     */
    {
        using Method = Experimental::MultiStage::RungeKuttaExplicitFourthOrder<Scalar>;
        auto ode = std::make_shared<ExponentialRhsODE>();
        auto method = std::make_shared<Method>();
        Experimental::ODESolver<ExponentialRhsODE> solverODE(ode, method);

        // initializes \f$ u=0 \f$ at the default independent-variable coordinate \f$ \xi=0 \f$
        Experimental::ODEVariables<Scalar> vars(0.0);
        solverODE.solve(vars, 1.0, 0.01);

        using std::exp;
        expectNear(vars.dofs(), exp(1.0) - 1.0, 1e-10, "Explicit ODE solve with rhs-only system failed");
    }

    /*
     * Solves the ODE: \f$ \frac{\mathrm{d}u(\xi)}{\mathrm{d}\xi} = -u(\xi) + \mathrm{e}^\xi \f$
     * on \f$ \xi \in [0,1] \f$ with \f$ u(0)=0 \f$
     * and a step size of \f$ \Delta \xi = 0.01 \f$
     *
     * note: `LinearScalarODE` must define `rhsJacobian()` because this test uses an implicit integration method.
     */
    {
        using Method = Experimental::MultiStage::Theta<Scalar>;
        auto ode = std::make_shared<LinearScalarODE>(-1.0);
        auto method = std::make_shared<Method>(0.5);
        Experimental::ODESolver<LinearScalarODE> solverODE(ode, method);

        // initializes \f$ u=0 \f$ at the default independent-variable coordinate \f$ \xi=0 \f$
        Experimental::ODEVariables<Scalar> vars(0.0);
        solverODE.solve(vars, 1.0, 0.01);

        using std::sinh;
        expectNear(vars.dofs(), sinh(1.0), 1e-5, "Implicit scalar ODE solve failed");
    }

    /*
     * Solves the ODE system:
     * \f[
     * \begin{aligned}
     * \frac{\mathrm{d}u_1(\xi)}{\mathrm{d}\xi} &= u_2(\xi), \\
     * \frac{\mathrm{d}u_2(\xi)}{\mathrm{d}\xi} &= -u_1(\xi)
     * \end{aligned}
     * \f]
     * with initial conditions \f$ u_1(0) = 1, \qquad u_2(0) = 0 \f$, \f$ \xi \in [0,0.5] \f$ and a step size of \f$ \Delta \xi = 0.01 \f$
     */
    {
        using Method = Experimental::MultiStage::Theta<Scalar>;
        using SolutionVector = HarmonicOscillatorODE::SolutionVector;

        auto ode = std::make_shared<HarmonicOscillatorODE>();
        auto method = std::make_shared<Method>(0.5);
        Experimental::ODESolver<HarmonicOscillatorODE> solverODE(ode, method);

        SolutionVector initial;
        initial[0] = 1.0;
        initial[1] = 0.0;
        // initializes \f$ u=0 \f$ at the default independent-variable coordinate \f$ \xi=0 \f$
        Experimental::ODEVariables<SolutionVector> vars(initial);
        solverODE.solve(vars, 0.5, 0.01);

        using std::cos;
        using std::sin;
        expectNear(vars.dofs()[0], cos(0.5), 1e-5, "Implicit vector ODE solve failed for the first component");
        expectNear(vars.dofs()[1], -sin(0.5), 1e-5, "Implicit vector ODE solve failed for the second component");
    }

    /*
     * Solves the ODE: \f$ \frac{\mathrm{d} (2u(\xi))}{\mathrm{d}\xi} = 1.0 \f$
     * on \f$ \xi \in [0,1] \f$ with \f$ u(0)=0 \f$
     * and a step size of \f$ \Delta \xi = 0.01 \f$
     *
     * note: by defining `derivativeQuantity` and `derivativeQuantityJacobian` in `ScaledDerivativeQuantityODE`, one can specify a custom derivative format, otherwise the form \f$ \frac{\mathrm{d}u(\xi)}{\mathrm{d}\xi} \f$ is implicitly assumed.
     */
    {
        using Method = Experimental::MultiStage::Theta<Scalar>;
        auto ode = std::make_shared<ScaledDerivativeQuantityODE>();
        auto method = std::make_shared<Method>(0.5);
        Experimental::ODESolver<ScaledDerivativeQuantityODE> solverODE(ode, method);

        // initializes \f$ u=0 \f$ at the default independent-variable coordinate \f$ \xi=0 \f$
        Experimental::ODEVariables<Scalar> vars(0.0);
        solverODE.solve(vars, 1.0, 0.01);

        expectNear(vars.dofs(), 0.5, 1e-10, "ODE solve with a custom differentiated quantity failed");
    }

    /*
     * Solves the ODE: \f$ \frac{\mathrm{d} (2u(\xi))}{\mathrm{d}\xi} = 1.0 \f$
     * on \f$ \xi \in [2,4] \f$ with \f$ u(2)=7 \f$
     * and a step size of \f$ \Delta \xi = 0.01 \f$
     *
     * note: here, a non-zero start for the integration is assigned
     */
    {
        using Method = Experimental::MultiStage::Theta<Scalar>;
        auto ode = std::make_shared<ScaledDerivativeQuantityODE>();
        auto method = std::make_shared<Method>(0.5);
        Experimental::ODESolver<ScaledDerivativeQuantityODE> solverODE(ode, method);

        using Level = Experimental::IndependentVariableLevel<Scalar>;
        // initializes vars at \f$ \xi = 2.0 \f$ with value 7.0
        Experimental::ODEVariables<Scalar> vars(7.0, Level{2.0});
        solverODE.solve(vars, 4.0, 0.01);

        expectNear(vars.dofs(), 8.0, 1e-10, "ODE solve with a custom initialization interval failed");
    }

    std::cout << "ODE solver tests passed" << std::endl;

    return 0;
}
