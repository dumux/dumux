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

class ExponentialRhsODE
{
public:
    using Scalar = double;
    using SolutionVector = Scalar;
    using ResidualType = Scalar;
    using JacobianMatrix = Scalar;
    using Variables = Experimental::Variables<SolutionVector>;

    void rhs(const Variables& vars, ResidualType& rhs) const
    {
        using std::exp;
        rhs = exp(vars.timeLevel().current());
    }
};

class LinearScalarODE
{
public:
    using Scalar = double;
    using SolutionVector = Scalar;
    using ResidualType = Scalar;
    using JacobianMatrix = Scalar;
    using Variables = Experimental::Variables<SolutionVector>;

    explicit LinearScalarODE(const Scalar lambda)
    : lambda_(lambda)
    {}

    void rhs(const Variables& vars, ResidualType& rhs) const
    {
        using std::exp;
        rhs = lambda_*vars.dofs() + exp(vars.timeLevel().current());
    }

    void jacobian(const Variables&, JacobianMatrix& jacobian) const
    { jacobian = lambda_; }

private:
    Scalar lambda_;
};

class HarmonicOscillatorODE
{
public:
    using Scalar = double;
    using SolutionVector = Dune::FieldVector<Scalar, 2>;
    using ResidualType = SolutionVector;
    using JacobianMatrix = Dune::FieldMatrix<Scalar, 2, 2>;
    using Variables = Experimental::Variables<SolutionVector>;

    void rhs(const Variables& vars, ResidualType& rhs) const
    {
        rhs[0] = vars.dofs()[1];
        rhs[1] = -vars.dofs()[0];
    }

    void jacobian(const Variables&, JacobianMatrix& jacobian) const
    {
        jacobian = 0.0;
        jacobian[0][1] = 1.0;
        jacobian[1][0] = -1.0;
    }
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

    {
        using Method = Experimental::MultiStage::RungeKuttaExplicitFourthOrder<Scalar>;
        auto ode = std::make_shared<ExponentialRhsODE>();
        auto method = std::make_shared<Method>();
        Experimental::ODESolver<ExponentialRhsODE> solver(ode, method);

        Experimental::Variables<Scalar> vars(0.0);
        solver.solve(vars, 1.0, 0.01);

        using std::exp;
        expectNear(vars.dofs(), exp(1.0) - 1.0, 1e-10,
                   "Explicit ODE solve with rhs-only system failed");
    }

    {
        using Method = Experimental::MultiStage::Theta<Scalar>;
        auto ode = std::make_shared<LinearScalarODE>(-1.0);
        auto method = std::make_shared<Method>(0.5);
        Experimental::ODESolver<LinearScalarODE> solver(ode, method);

        Experimental::Variables<Scalar> vars(0.0);
        solver.solve(vars, 1.0, 0.01);

        using std::sinh;
        expectNear(vars.dofs(), sinh(1.0), 1e-5,
                   "Implicit scalar ODE solve failed");
    }

    {
        using Method = Experimental::MultiStage::Theta<Scalar>;
        using SolutionVector = HarmonicOscillatorODE::SolutionVector;

        auto ode = std::make_shared<HarmonicOscillatorODE>();
        auto method = std::make_shared<Method>(0.5);
        Experimental::ODESolver<HarmonicOscillatorODE> solver(ode, method);

        SolutionVector initial;
        initial[0] = 1.0;
        initial[1] = 0.0;
        Experimental::Variables<SolutionVector> vars(initial);
        solver.solve(vars, 0.5, 0.01);

        using std::cos;
        using std::sin;
        expectNear(vars.dofs()[0], cos(0.5), 1e-5,
                   "Implicit vector ODE solve failed for the first component");
        expectNear(vars.dofs()[1], -sin(0.5), 1e-5,
                   "Implicit vector ODE solve failed for the second component");
    }

    std::cout << "ODE solver tests passed" << std::endl;

    return 0;
}
