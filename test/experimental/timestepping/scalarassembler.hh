// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Experimental
 * \brief Assembler and linear solver for the scalar test equation du/dt = lambda*u + exp(t)
 *        of the time stepping tests
 */
#ifndef DUMUX_TEST_TIMESTEPPING_SCALAR_ASSEMBLER_HH
#define DUMUX_TEST_TIMESTEPPING_SCALAR_ASSEMBLER_HH

#include <cmath>
#include <optional>

#include <dune/common/exceptions.hh>

#include <dumux/experimental/common/variables.hh>
#include <dumux/experimental/timestepping/timelevel.hh>
#include <dumux/experimental/timestepping/implicitsolve.hh>

namespace Dumux {

class ScalarAssembler
{
public:
    using Scalar = double;
    using SolutionVector = Scalar;
    using ResidualType = Scalar;
    using JacobianMatrix = Scalar;
    using Variables = Experimental::Variables<SolutionVector>;
    using SolveParams = Experimental::ImplicitSolveParams<Scalar, ResidualType>;

    explicit ScalarAssembler(Scalar lambda = 0.0) : lambda_(lambda) {}

    void setLinearSystem() {}
    JacobianMatrix& jacobian() { return jac_; }
    ResidualType& residual() { return res_; }

    void assembleResidual(const Variables& vars)
    {
        if (!params_)
            DUNE_THROW(Dune::InvalidStateException, "No implicit solve set");

        const auto& weights = params_->weights;
        res_ = weights.storage*storage_(vars) + weights.spatial*spatial_(vars);
        if (params_->constant)
            res_ += *params_->constant;
    }

    void assembleJacobianAndResidual(const Variables& vars)
    {
        assembleResidual(vars);
        jac_ = params_->weights.storage - lambda_*params_->weights.spatial;
    }

    void resetTimeStep(const SolutionVector&) {}

    void setImplicitSolve(Variables& vars, const SolveParams& params)
    {
        vars.updateTime(params.timeLevel);
        params_ = params;
    }

    void assembleOperators(Variables& vars,
                           const Experimental::TimeLevel<Scalar>& timeLevel,
                           Experimental::OperatorEvaluations<ResidualType> evaluations)
    {
        vars.updateTime(timeLevel);
        if (evaluations.storage)
            *evaluations.storage = storage_(vars);
        if (evaluations.spatial)
            *evaluations.spatial = spatial_(vars);
    }

private:
    Scalar storage_(const Variables& vars) const
    { return vars.dofs(); }

    Scalar spatial_(const Variables& vars) const
    {
        using std::exp;
        return -(lambda_*vars.dofs() + exp(vars.timeLevel().current()));
    }

    Scalar lambda_;
    ResidualType res_;
    JacobianMatrix jac_;
    std::optional<SolveParams> params_;
};

class ScalarLinearSolver
{
public:
    void setResidualReduction(double residualReduction) {}

    bool solve(const double& A, double& x, const double& b) const
    { x = b/A; return true; }

    double norm(const double residual) const
    {
        using std::abs;
        return abs(residual);
    }
};

} // end namespace Dumux

#endif
