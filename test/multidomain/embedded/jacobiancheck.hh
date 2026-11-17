// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief Helpers to check the Jacobian of a coupled multidomain system
 */
#ifndef DUMUX_TEST_MULTIDOMAIN_EMBEDDED_JACOBIANCHECK_HH
#define DUMUX_TEST_MULTIDOMAIN_EMBEDDED_JACOBIANCHECK_HH

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <random>
#include <utility>
#include <vector>

#include <dune/common/hybridutilities.hh>
#include <dune/common/indices.hh>

namespace Dumux::Test {

//! Call f(domainIdx, dofIdx, priVarIdx, value) for every entry of a multidomain vector
template<class MultiDomainVector, class F>
void forEachEntry(MultiDomainVector& v, F&& f)
{
    using namespace Dune::Hybrid;
    forEach(integralRange(Dune::Hybrid::size(v)), [&](const auto domainIdx)
    {
        auto& subVector = v[domainIdx];
        for (std::size_t dofIdx = 0; dofIdx < subVector.size(); ++dofIdx)
            for (std::size_t pvIdx = 0; pvIdx < subVector[dofIdx].size(); ++pvIdx)
                f(std::size_t(domainIdx), dofIdx, pvIdx, subVector[dofIdx][pvIdx]);
    });
}

//! The Euclidean norm of a multidomain vector
template<class MultiDomainVector>
double norm(const MultiDomainVector& v)
{
    double sum = 0.0;
    forEachEntry(v, [&](std::size_t, std::size_t, std::size_t, const double x) { sum += x*x; });
    return std::sqrt(sum);
}

//! The result of the comparison of the product of the Jacobian with a vector and a finite difference of the residual
struct JacobianCheckEntry
{
    double epsilon; //!< the size of the perturbation of the solution
    double relativeError; //!< norm of the difference relative to the norm of J v
};

/*!
 * \brief Compare J v to (R(u + eps v) - R(u))/eps for a random direction v
 *
 * \param assembler the multidomain assembler
 * \param couplingManager the coupling manager of the assembler
 * \param sol the solution u at which the Jacobian is evaluated
 * \param epsilons the sizes of the perturbation
 * \param isConstrained constrained(domainIdx, dofIdx) is true for dofs where the direction is zero
 * \param seed the seed of the random number generator for the direction
 *
 * The direction is a unit vector with random components. The assembler and the
 * coupling manager are in the state of the solution u when the function returns.
 * If an entry of the Jacobian is missing or has the wrong sign, the relative error
 * does not go to zero for any epsilon. If the residual is linear in the solution,
 * the relative error is the round-off error of the residual divided by epsilon.
 */
template<class Assembler, class CouplingManager, class SolutionVector>
std::vector<JacobianCheckEntry>
compareJacobianWithFiniteDifferences(Assembler& assembler,
                                     CouplingManager& couplingManager,
                                     const SolutionVector& sol,
                                     const std::vector<double>& epsilons,
                                     const std::function<bool(std::size_t, std::size_t)>& isConstrained = {},
                                     unsigned int seed = 1)
{
    couplingManager.updateSolution(sol);
    assembler.assembleJacobianAndResidual(sol);
    const auto residual0 = assembler.residual();

    auto direction = sol;
    std::mt19937 generator(seed);
    std::uniform_real_distribution<double> distribution(-1.0, 1.0);
    forEachEntry(direction, [&](std::size_t domainIdx, std::size_t dofIdx, std::size_t, double& x)
    { x = (isConstrained && isConstrained(domainIdx, dofIdx)) ? 0.0 : distribution(generator); });
    const auto directionNorm = norm(direction);
    forEachEntry(direction, [&](std::size_t, std::size_t, std::size_t, double& x) { x /= directionNorm; });

    auto jacobianTimesDirection = residual0;
    jacobianTimesDirection = 0.0;
    assembler.jacobian().mv(direction, jacobianTimesDirection);
    const auto referenceNorm = norm(jacobianTimesDirection);

    std::vector<JacobianCheckEntry> result;
    for (const auto epsilon : epsilons)
    {
        auto perturbedSol = sol;
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(perturbedSol)), [&](const auto domainIdx)
        {
            for (std::size_t dofIdx = 0; dofIdx < perturbedSol[domainIdx].size(); ++dofIdx)
                for (std::size_t pvIdx = 0; pvIdx < perturbedSol[domainIdx][dofIdx].size(); ++pvIdx)
                    perturbedSol[domainIdx][dofIdx][pvIdx] += epsilon*direction[domainIdx][dofIdx][pvIdx];
        });

        couplingManager.updateSolution(perturbedSol);
        assembler.assembleResidual(perturbedSol);
        auto difference = assembler.residual();
        difference -= residual0;
        difference /= epsilon;
        difference -= jacobianTimesDirection;
        result.push_back({epsilon, norm(difference)/referenceNorm});
    }

    couplingManager.updateSolution(sol);
    assembler.updateGridVariables(sol);
    return result;
}

//! The largest absolute difference of the entries of two multidomain matrices with the same pattern and the largest absolute entry
struct MatrixDifference
{
    double maxAbsDifference = 0.0;
    double maxAbsEntry = 0.0;
    bool identical = true; //!< the entries are equal bit by bit
};

template<class MultiDomainMatrix>
MatrixDifference compareMatrices(const MultiDomainMatrix& a, const MultiDomainMatrix& b)
{
    MatrixDifference result;
    using namespace Dune::Hybrid;
    forEach(integralRange(Dune::index_constant<MultiDomainMatrix::N()>{}), [&](const auto i)
    {
        forEach(integralRange(Dune::index_constant<MultiDomainMatrix::M()>{}), [&](const auto j)
        {
            const auto& blockMatrixA = a[i][j];
            const auto& blockMatrixB = b[i][j];
            for (auto row = blockMatrixA.begin(); row != blockMatrixA.end(); ++row)
            {
                for (auto col = row->begin(); col != row->end(); ++col)
                {
                    if (!blockMatrixB.exists(row.index(), col.index()))
                    {
                        result.identical = false;
                        result.maxAbsDifference = std::numeric_limits<double>::infinity();
                        continue;
                    }
                    const auto& blockA = *col;
                    const auto& blockB = blockMatrixB[row.index()][col.index()];
                    for (std::size_t r = 0; r < blockA.N(); ++r)
                    {
                        for (std::size_t c = 0; c < blockA.M(); ++c)
                        {
                            result.maxAbsEntry = std::max(result.maxAbsEntry, std::abs(blockA[r][c]));
                            result.maxAbsDifference = std::max(result.maxAbsDifference, std::abs(blockA[r][c] - blockB[r][c]));
                            if (blockA[r][c] != blockB[r][c])
                                result.identical = false;
                        }
                    }
                }
            }
        });
    });
    return result;
}

} // end namespace Dumux::Test

#endif
