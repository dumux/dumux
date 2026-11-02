// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Compare the weighted assembly of the multi-stage assembler with the assembly
 *        of the instationary finite volume assembler for implicit Euler
 *
 * With the weights \f$ w_M = 1 \f$, \f$ w_R = \Delta t \f$ and the constant
 * \f$ c = -M(x_n) \f$, the multi-stage residual \f$ M(x) + \Delta t R(x) - M(x_n) \f$ is
 * \f$ \Delta t \f$ times the residual of the instationary assembler,
 * \f$ (M(x) - M(x_n))/\Delta t + R(x) \f$, and so is its Jacobian.
 */
#ifndef DUMUX_TEST_ASSEMBLY_WEIGHTED_ASSEMBLY_HH
#define DUMUX_TEST_ASSEMBLY_WEIGHTED_ASSEMBLY_HH

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>

#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/multistageassembler.hh>
#include <dumux/experimental/timestepping/multistagemethods.hh>

#include <test/assembly/constraintrows.hh>

namespace Dumux::Test {

struct WeightedAssemblyDifferences
{
    double residual = 0.0; //!< largest difference of the residuals relative to their largest entry
    double jacobian = 0.0; //!< largest difference of the Jacobians relative to their largest entry
    double residualOnly = 0.0; //!< largest difference of the residuals assembled without and with the Jacobian
    double constraint = 0.0; //!< largest difference in the constraint rows, which are not scaled by the time step size
    std::size_t numCompared = 0; //!< number of compared rows
    std::size_t numConstraintRows = 0; //!< number of compared constraint rows
    std::size_t numPeriodicRows = 0; //!< number of compared constraint rows of periodic degrees of freedom
};

/*!
 * \brief Assemble one implicit Euler step with both assemblers
 * \param x The state at which the step is assembled
 * \param xOld The state at the beginning of the step
 */
template<class TypeTag, class Problem, class GridGeometry, class GridVariables, class SolutionVector>
WeightedAssemblyDifferences compareWeightedAssembly(std::shared_ptr<Problem> problem,
                                                    std::shared_ptr<GridGeometry> gridGeometry,
                                                    std::shared_ptr<GridVariables> gridVariables,
                                                    const SolutionVector& x,
                                                    const SolutionVector& xOld,
                                                    double dt)
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using std::abs; using std::max;

    // the instationary finite volume assembler
    gridVariables->init(xOld);
    gridVariables->update(x);
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, 10.0*dt);
    using FVAssembler_ = FVAssembler<TypeTag, DiffMethod::numeric>;
    FVAssembler_ fvAssembler(problem, gridGeometry, gridVariables, timeLoop, xOld);
    fvAssembler.setLinearSystem();
    fvAssembler.assembleJacobianAndResidual(x);

    // the multi-stage assembler
    using MultiStageAssembler_ = Experimental::MultiStageAssembler<TypeTag, DiffMethod::numeric>;
    auto method = std::make_shared<Experimental::MultiStage::ImplicitEuler<Scalar>>();
    MultiStageAssembler_ msAssembler(problem, gridGeometry, gridVariables, method, xOld);
    msAssembler.setLinearSystem();

    using ResidualType = typename MultiStageAssembler_::ResidualType;
    auto storageAtOld = std::make_shared<ResidualType>();
    gridVariables->update(xOld);
    msAssembler.assembleOperators(xOld, Experimental::TimeLevel<Scalar>(0.0), {.storage = storageAtOld.get()});
    *storageAtOld *= -1.0;
    gridVariables->update(x);

    typename MultiStageAssembler_::SolveParams params{
        Experimental::TimeLevel<Scalar>(dt, 0.0, 1.0), {.inertia = 0.0, .storage = 1.0, .spatial = dt}, storageAtOld
    };
    msAssembler.setImplicitSolve(x, params);
    msAssembler.assembleJacobianAndResidual(x);

    WeightedAssemblyDifferences result;

    const auto residualWithJacobian = msAssembler.residual();
    msAssembler.assembleResidual(x);
    auto residualDifference = residualWithJacobian;
    residualDifference -= msAssembler.residual();
    result.residualOnly = residualDifference.infinity_norm()/max(residualWithJacobian.infinity_norm(), 1e-300);

    // rows with a constraint are not scaled by the time step size
    const auto& jacOld = fvAssembler.jacobian();
    const auto& jacMs = msAssembler.jacobian();
    const auto isConstraintRow = [&] (std::size_t i) { return Test::isConstraintRow(jacOld[i], i); };
    const auto isPeriodicRow = [&] (std::size_t i) { return Test::isConstraintRow(jacOld[i], i, /*periodic=*/true); };

    double residualScale = 0.0, jacobianScale = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        if (isConstraintRow(i))
            continue;
        for (int r = 0; r < int(x[i].size()); ++r)
            residualScale = max(residualScale, abs(dt*fvAssembler.residual()[i][r]));
        const auto& row = jacOld[i];
        for (auto col = row.begin(); col != row.end(); ++col)
            for (int r = 0; r < int(col->N()); ++r)
                for (int c = 0; c < int(col->M()); ++c)
                    jacobianScale = max(jacobianScale, abs(dt*(*col)[r][c]));
    }

    for (std::size_t i = 0; i < x.size(); ++i)
    {
        if (isConstraintRow(i))
        {
            ++result.numConstraintRows;
            if (isPeriodicRow(i))
                ++result.numPeriodicRows;
            for (int r = 0; r < int(x[i].size()); ++r)
                result.constraint = max(result.constraint, abs(msAssembler.residual()[i][r] - fvAssembler.residual()[i][r]));
            const auto& row = jacOld[i];
            for (auto col = row.begin(); col != row.end(); ++col)
                for (int r = 0; r < int(col->N()); ++r)
                    for (int c = 0; c < int(col->M()); ++c)
                        result.constraint = max(result.constraint, abs(jacMs[i][col.index()][r][c] - (*col)[r][c]));
            continue;
        }

        ++result.numCompared;
        for (int r = 0; r < int(x[i].size()); ++r)
            result.residual = max(result.residual,
                                  abs(msAssembler.residual()[i][r] - dt*fvAssembler.residual()[i][r])/residualScale);

        const auto& row = jacOld[i];
        for (auto col = row.begin(); col != row.end(); ++col)
            for (int r = 0; r < int(col->N()); ++r)
                for (int c = 0; c < int(col->M()); ++c)
                    result.jacobian = max(result.jacobian,
                                          abs(jacMs[i][col.index()][r][c] - dt*(*col)[r][c])/jacobianScale);
    }

    return result;
}

} // end namespace Dumux::Test

#endif
