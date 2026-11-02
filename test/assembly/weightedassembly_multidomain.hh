// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Compare the weighted assembly of the multi-stage multi-domain assembler with the assembly
 *        of the instationary multi-domain assembler for implicit Euler, including the coupling blocks
 *
 * With the weights \f$ w_M = 1 \f$, \f$ w_R = \Delta t \f$ and the constant
 * \f$ c = -M(x_n) \f$, the multi-stage residual \f$ M(x) + \Delta t R(x) - M(x_n) \f$ is
 * \f$ \Delta t \f$ times the residual of the instationary assembler,
 * \f$ (M(x) - M(x_n))/\Delta t + R(x) \f$, and so is its Jacobian.
 * Constraint rows are not scaled by the time step size.
 */
#ifndef DUMUX_TEST_ASSEMBLY_WEIGHTED_ASSEMBLY_MULTIDOMAIN_HH
#define DUMUX_TEST_ASSEMBLY_WEIGHTED_ASSEMBLY_MULTIDOMAIN_HH

#include <algorithm>
#include <cmath>
#include <memory>

#include <dune/common/hybridutilities.hh>

#include <dumux/experimental/timestepping/timelevel.hh>

#include <test/assembly/constraintrows.hh>

namespace Dumux::Test {

struct MultiDomainWeightedAssemblyDifferences
{
    double residual = 0.0; //!< largest scaled difference of the residuals in the unconstrained rows
    double jacobian = 0.0; //!< largest scaled difference of the Jacobian blocks in the unconstrained rows
    double coupling = 0.0; //!< largest scaled difference of the off-diagonal Jacobian blocks in the unconstrained rows
    double constraint = 0.0; //!< largest difference in the constraint rows
    std::size_t numCompared = 0; //!< number of compared rows
    std::size_t numCouplingEntries = 0; //!< number of compared nonzero entries of off-diagonal blocks
    std::size_t numConstraintRows = 0; //!< number of compared constraint rows
    std::size_t numPeriodicRows = 0; //!< number of compared constraint rows of periodic degrees of freedom
};

/*!
 * \brief Assemble one implicit Euler step with both assemblers
 * \param oldAssembler The instationary assembler (previous solution xOld)
 * \param msAssembler The multi-stage assembler
 * \param x The state at which the step is assembled
 * \param dt The time step size
 * \note The grid variables are expected to be initialized with the solution at the beginning of the step
 */
template<class OldAssembler, class MultiStageAssembler, class SolutionVector>
MultiDomainWeightedAssemblyDifferences compareWeightedMultiDomainAssembly(OldAssembler& oldAssembler,
                                                                          MultiStageAssembler& msAssembler,
                                                                          const SolutionVector& x,
                                                                          const SolutionVector& xOld,
                                                                          double dt)
{
    using std::abs; using std::max;
    using Scalar = typename MultiStageAssembler::Scalar;

    oldAssembler.setLinearSystem();
    oldAssembler.updateGridVariables(x);
    oldAssembler.assembleJacobianAndResidual(x);

    msAssembler.setLinearSystem();
    using ResidualType = typename MultiStageAssembler::ResidualType;
    auto storageAtOld = std::make_shared<ResidualType>();
    msAssembler.updateGridVariables(xOld);
    msAssembler.assembleOperators(xOld, Experimental::TimeLevel<Scalar>(0.0), {.storage = storageAtOld.get()});
    *storageAtOld *= -1.0;
    msAssembler.updateGridVariables(x);

    typename MultiStageAssembler::SolveParams params{
        Experimental::TimeLevel<Scalar>(dt, 0.0, 1.0), {.inertia = 0.0, .storage = 1.0, .spatial = dt}, storageAtOld
    };
    msAssembler.setImplicitSolve(x, params);
    msAssembler.assembleJacobianAndResidual(x);

    const auto& jacOld = oldAssembler.jacobian();
    const auto& jacMs = msAssembler.jacobian();
    const auto& resOld = oldAssembler.residual();
    const auto& resMs = msAssembler.residual();

    // a constraint row is a constraint row of its diagonal block and zero in all other blocks
    const auto isConstraintRow = [&] (auto domainI, std::size_t i, bool periodic = false)
    {
        bool result = Test::isConstraintRow(jacOld[domainI][domainI][i], i, periodic);
        Dune::Hybrid::forEach(std::make_index_sequence<std::decay_t<decltype(jacOld)>::N()>{}, [&](auto domainJ)
        {
            if constexpr (decltype(domainJ)::value != decltype(domainI)::value)
            {
                const auto& row = jacOld[domainI][domainJ][i];
                for (auto col = row.begin(); col != row.end(); ++col)
                    for (int r = 0; r < int(col->N()); ++r)
                        for (int c = 0; c < int(col->M()); ++c)
                            if ((*col)[r][c] != 0.0)
                                result = false;
            }
        });
        return result;
    };

    MultiDomainWeightedAssemblyDifferences result;
    double residualScale = 0.0, jacobianScale = 0.0;
    constexpr auto numDomains = std::decay_t<decltype(jacOld)>::N();

    Dune::Hybrid::forEach(std::make_index_sequence<numDomains>{}, [&](auto domainI)
    {
        for (std::size_t i = 0; i < resOld[domainI].size(); ++i)
        {
            if (isConstraintRow(domainI, i))
                continue;
            for (int r = 0; r < int(resOld[domainI][i].size()); ++r)
                residualScale = max(residualScale, abs(dt*resOld[domainI][i][r]));
            Dune::Hybrid::forEach(std::make_index_sequence<numDomains>{}, [&](auto domainJ)
            {
                const auto& row = jacOld[domainI][domainJ][i];
                for (auto col = row.begin(); col != row.end(); ++col)
                    for (int r = 0; r < int(col->N()); ++r)
                        for (int c = 0; c < int(col->M()); ++c)
                            jacobianScale = max(jacobianScale, abs(dt*(*col)[r][c]));
            });
        }
    });

    Dune::Hybrid::forEach(std::make_index_sequence<numDomains>{}, [&](auto domainI)
    {
        for (std::size_t i = 0; i < resOld[domainI].size(); ++i)
        {
            const bool constraint = isConstraintRow(domainI, i);
            if (constraint)
            {
                ++result.numConstraintRows;
                if (isConstraintRow(domainI, i, /*periodic=*/true))
                    ++result.numPeriodicRows;
            }
            else
                ++result.numCompared;

            for (int r = 0; r < int(resOld[domainI][i].size()); ++r)
            {
                if (constraint)
                    result.constraint = max(result.constraint, abs(resMs[domainI][i][r] - resOld[domainI][i][r]));
                else
                    result.residual = max(result.residual, abs(resMs[domainI][i][r] - dt*resOld[domainI][i][r])/residualScale);
            }

            Dune::Hybrid::forEach(std::make_index_sequence<numDomains>{}, [&](auto domainJ)
            {
                const auto& row = jacOld[domainI][domainJ][i];
                for (auto col = row.begin(); col != row.end(); ++col)
                    for (int r = 0; r < int(col->N()); ++r)
                        for (int c = 0; c < int(col->M()); ++c)
                        {
                            const auto ms = jacMs[domainI][domainJ][i][col.index()][r][c];
                            if (constraint)
                                result.constraint = max(result.constraint, abs(ms - (*col)[r][c]));
                            else
                            {
                                const auto difference = abs(ms - dt*(*col)[r][c])/jacobianScale;
                                result.jacobian = max(result.jacobian, difference);
                                if (domainI != domainJ)
                                {
                                    result.coupling = max(result.coupling, difference);
                                    if ((*col)[r][c] != 0.0)
                                        ++result.numCouplingEntries;
                                }
                            }
                        }
            });
        }
    });

    return result;
}

} // end namespace Dumux::Test

#endif
