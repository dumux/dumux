// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#ifndef DUMUX_TEST_NEWTON_COMMON_HH
#define DUMUX_TEST_NEWTON_COMMON_HH

#include <cmath>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

/*
  This corresponding test currently solves a scalar non-linear equation (
  or small system of nonlinear equations) using Newton's method.
  The Mock classes expose which dependencies the
  current implementation has on different other classes it interacts with.
  This test is to ensure that the dependencies do not grow more in the future.
  The Mock classes can be step-wise reduced in complexity once the Newton
  implementation required a smaller interface from assembler and solver.
 */

namespace Dumux {

// R(x) = x^2 - 5 = 0
class MockScalarAssembler
{
public:
    using Scalar = double;
    using ResidualType = Scalar;
    using JacobianMatrix = Scalar;
    using SolutionVector = Scalar;
    using Variables = Scalar;

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

private:
    JacobianMatrix jac_;
    ResidualType res_;
};

class MockScalarLinearSolver
{
public:
    void setResidualReduction(double residualReduction) {}

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

class MockSystemValuedAssembler
{
public:
    using Scalar = double;
    using ResidualType = Dune::FieldVector<Scalar, 2>;
    using JacobianMatrix = Dune::FieldMatrix<Scalar, 2, 2>;
    using SolutionVector = Dune::FieldVector<Scalar, 2>;
    using Variables = Dune::FieldVector<Scalar, 2>;

    void setLinearSystem() {}

    void assembleResidual(const ResidualType& sol)
    {
        for (int i = 0; i < 2; ++i)
            res_[i] = sol[i]*sol[i] - 5.0;
    }

    void assembleJacobianAndResidual (const ResidualType& sol)
    {
        assembleResidual(sol);
        jac_ = 2.0*sol;
    }

    JacobianMatrix& jacobian() { return jac_; }

    ResidualType& residual() { return res_; }

private:
    JacobianMatrix jac_;
    ResidualType res_;
};

class MockSystemLinearSolver
{
public:
    void setResidualReduction(double residualReduction) {}

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

} // end namespace Dumux

#endif
