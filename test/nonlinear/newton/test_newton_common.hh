// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#ifndef DUMUX_TEST_NEWTON_COMMON_HH
#define DUMUX_TEST_NEWTON_COMMON_HH

#include <cmath>
#include <functional>
#include <utility>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

/*
  Mock assemblers and linear solvers for solving a scalar nonlinear equation
  or a small system of nonlinear equations with Newton's method.
  The Mock classes expose which dependencies the
  current implementation has on different other classes it interacts with.
  The tests ensure that the dependencies do not grow more in the future.
  The Mock classes can be step-wise reduced in complexity once the Newton
  implementation required a smaller interface from assembler and solver.
 */

namespace Dumux {

// R(x) = f(x) = 0, by default with f(x) = x^2 - 5
class MockScalarAssembler
{
public:
    using Scalar = double;
    using ResidualType = Scalar;
    using JacobianMatrix = Scalar;
    using SolutionVector = Scalar;
    using Variables = Scalar;

    MockScalarAssembler()
    : MockScalarAssembler([](Scalar x){ return x*x - 5.0; }, [](Scalar x){ return 2.0*x; })
    {}

    MockScalarAssembler(std::function<Scalar(Scalar)> f, std::function<Scalar(Scalar)> df)
    : f_(std::move(f)), df_(std::move(df))
    {}

    void setLinearSystem() {}

    void assembleResidual(const ResidualType& sol)
    {
        res_ = f_(sol);
    }

    void assembleJacobianAndResidual (const ResidualType& sol)
    {
        assembleResidual(sol);
        jac_ = df_(sol);
    }

    JacobianMatrix& jacobian() { return jac_; }

    ResidualType& residual() { return res_; }

private:
    std::function<Scalar(Scalar)> f_;
    std::function<Scalar(Scalar)> df_;
    JacobianMatrix jac_;
    ResidualType res_;
};

class MockScalarLinearSolver
{
public:
    void setResidualReduction(double) {}

    bool solve(const double& A, double& x, const double& b) const
    {
        x = b/A;
        return true;
    }

    double norm(const double& residual) const
    {
        using std::abs;
        return abs(residual);
    }
};

// R_i(x) = f(x_i) = 0, i = 0, 1
class MockSystemAssembler
{
public:
    using Scalar = double;
    using ResidualType = Dune::FieldVector<Scalar, 2>;
    using JacobianMatrix = Dune::FieldMatrix<Scalar, 2, 2>;
    using SolutionVector = Dune::FieldVector<Scalar, 2>;
    using Variables = Dune::FieldVector<Scalar, 2>;

    MockSystemAssembler(std::function<Scalar(Scalar)> f, std::function<Scalar(Scalar)> df)
    : f_(std::move(f)), df_(std::move(df))
    {}

    void setLinearSystem() {}

    void assembleResidual(const SolutionVector& sol)
    {
        for (int i = 0; i < 2; ++i)
            res_[i] = f_(sol[i]);
    }

    void assembleJacobianAndResidual (const SolutionVector& sol)
    {
        assembleResidual(sol);
        jac_ = 0.0;
        for (int i = 0; i < 2; ++i)
            jac_[i][i] = df_(sol[i]);
    }

    JacobianMatrix& jacobian() { return jac_; }

    ResidualType& residual() { return res_; }

private:
    std::function<Scalar(Scalar)> f_;
    std::function<Scalar(Scalar)> df_;
    JacobianMatrix jac_;
    ResidualType res_;
};

class MockSystemLinearSolver
{
public:
    using Matrix = Dune::FieldMatrix<double, 2, 2>;
    using Vector = Dune::FieldVector<double, 2>;

    void setResidualReduction(double) {}

    bool solve(const Matrix& A, Vector& x, const Vector& b) const
    {
        A.solve(x, b);
        return true;
    }

    double norm(const Vector& residual) const
    {
        return residual.two_norm();
    }
};

} // end namespace Dumux

#endif
