// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Experimental
 * \brief Check of the contract between time steppers and assemblers (see ImplicitSolveParams)
 */
#ifndef DUMUX_TEST_TIMESTEPPING_IMPLICIT_SOLVE_CHECK_HH
#define DUMUX_TEST_TIMESTEPPING_IMPLICIT_SOLVE_CHECK_HH

#include <cmath>
#include <memory>
#include <random>
#include <iostream>
#include <algorithm>

#include <dune/common/exceptions.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/typetraits.hh>

#include <dumux/experimental/timestepping/timelevel.hh>
#include <dumux/experimental/timestepping/implicitsolve.hh>

namespace Dumux::Test {

namespace Detail {

//! Call f(a, b, c, d) for all scalar entries of four vectors with the same block structure
template<class V, class F>
void forEachEntry(V& a, V& b, V& c, V& d, F&& f)
{
    if constexpr (Dune::IsNumber<std::decay_t<V>>::value)
        f(a, b, c, d);
    else if constexpr (requires { a.size(); a[0]; })
    {
        for (std::size_t i = 0; i < a.size(); ++i)
            forEachEntry(a[i], b[i], c[i], d[i], f);
    }
    else
        Dune::Hybrid::forEach(std::make_index_sequence<std::decay_t<V>::size()>{}, [&](auto i)
        { forEachEntry(a[i], b[i], c[i], d[i], f); });
}

} // end namespace Detail

/*!
 * \brief Checks that an assembler implements the contract of an implicit solve
 *
 * For a given solution, the residual assembled for the weights (wI, wM, wR) and a
 * constant c has to equal wI I + wM M + wR R + c in all unconstrained equations, where
 * I, M and R are the operators returned by assembleOperators, and has to be independent
 * of c in the constrained equations. The inertia operator has to vanish for problems
 * without inertia.
 * \return the number of constrained equations
 */
template<class Assembler, class SolutionVector, class Scalar>
std::size_t checkImplicitSolveContract(Assembler& assembler, const SolutionVector& x, const Scalar time,
                                       const Scalar storageWeight, const Scalar spatialWeight,
                                       const Scalar inertiaWeight = 0.0)
{
    using Residual = typename Assembler::ResidualType;
    using SolveParams = Experimental::ImplicitSolveParams<Scalar, Residual>;
    const Experimental::TimeLevel<Scalar> timeLevel(time);

    Residual inertia, storage, spatial;
    assembler.assembleOperators(x, timeLevel, {.inertia = &inertia, .storage = &storage, .spatial = &spatial});

    constexpr bool hasInertia = [] {
        if constexpr (requires { Assembler::hasInertia(); })
            return Assembler::hasInertia();
        else
            return false;
    }();
    if constexpr (!hasInertia)
    {
        Scalar maxInertia = 0.0;
        Detail::forEachEntry(inertia, inertia, inertia, inertia, [&](auto& i, auto&, auto&, auto&)
        { using std::abs; using std::max; maxInertia = max(maxInertia, Scalar(abs(i))); });
        if (maxInertia != 0.0)
            DUNE_THROW(Dune::Exception, "The inertia operator of a problem without inertia is not zero");
    }

    auto constant = std::make_shared<Residual>(storage);
    std::size_t counter = 0;
    Detail::forEachEntry(*constant, storage, storage, storage, [&](auto& c, auto&, auto&, auto&)
    { c = std::sin(1.0 + 0.37*(counter++)); });

    assembler.setImplicitSolve(x, SolveParams{timeLevel, {inertiaWeight, storageWeight, spatialWeight}, nullptr});
    assembler.assembleJacobianAndResidual(x);
    Residual withoutConstant = assembler.residual();

    assembler.setImplicitSolve(x, SolveParams{timeLevel, {inertiaWeight, storageWeight, spatialWeight}, constant});
    assembler.assembleJacobianAndResidual(x);
    Residual withConstant = assembler.residual();

    Scalar scale = 0.0;
    Detail::forEachEntry(withoutConstant, storage, spatial, inertia, [&](auto& r, auto& m, auto& s, auto& i)
    {
        using std::abs; using std::max;
        scale = max({scale, abs(r), abs(storageWeight*m), abs(spatialWeight*s), abs(inertiaWeight*i)});
    });

    std::size_t numConstrained = 0;
    Scalar maxConstantError = 0.0, maxOperatorError = 0.0;
    Detail::forEachEntry(withoutConstant, withConstant, *constant, storage, [&](auto& r0, auto& r, auto& c, auto&)
    {
        using std::abs; using std::max;
        if (r == r0)
            ++numConstrained;
        else
            maxConstantError = max(maxConstantError, abs((r - r0) - c));
    });

    Residual weighted = storage;
    weighted *= storageWeight;
    weighted.axpy(spatialWeight, spatial);
    weighted.axpy(inertiaWeight, inertia);
    Detail::forEachEntry(withoutConstant, withConstant, weighted, storage, [&](auto& r0, auto& r, auto& w, auto&)
    {
        using std::abs; using std::max;
        if (r != r0)
            maxOperatorError = max(maxOperatorError, abs(r0 - w));
    });

    std::cout << "[ImplicitSolveContract] constrained equations: " << numConstrained
              << ", largest deviation of the added constant: " << maxConstantError
              << ", largest deviation from the weighted operators: " << maxOperatorError/scale
              << " (relative to " << scale << ")" << std::endl;

    if (maxConstantError > 1e-12*std::max(scale, 1.0))
        DUNE_THROW(Dune::Exception, "The constant of the implicit solve is not added to the unconstrained equations");
    if (maxOperatorError > 1e-10*scale)
        DUNE_THROW(Dune::Exception, "The residual does not equal the weighted sum of the operator evaluations");

    return numConstrained;
}

/*!
 * \brief The largest deviation of the assembled Jacobian from central differences of the residual
 *        in a random direction, relative to the largest entry of the Jacobian-vector product
 * \note For single-domain assemblers with block vectors
 */
template<class Assembler, class SolutionVector, class Scalar>
Scalar jacobianDeviationFromFiniteDifferences(Assembler& assembler, const SolutionVector& x, const Scalar time,
                                              const Experimental::OperatorWeights<Scalar>& weights,
                                              const Scalar epsilon, const unsigned int seed = 1)
{
    using Residual = typename Assembler::ResidualType;
    using SolveParams = Experimental::ImplicitSolveParams<Scalar, Residual>;
    assembler.setImplicitSolve(x, SolveParams{Experimental::TimeLevel<Scalar>(time), weights, nullptr});

    const auto residualAt = [&](const SolutionVector& y)
    {
        assembler.updateGridVariables(y);
        assembler.assembleJacobianAndResidual(y);
        return Residual(assembler.residual());
    };

    const auto residual = residualAt(x);
    const auto jacobian = assembler.jacobian();

    std::mt19937 generator(seed);
    std::uniform_real_distribution<Scalar> distribution(-1.0, 1.0);
    Residual direction(residual);
    for (auto& block : direction)
        for (auto& entry : block)
            entry = distribution(generator);

    Residual jacobianTimesDirection(residual);
    jacobian.mv(direction, jacobianTimesDirection);

    auto xPlus = x, xMinus = x;
    for (std::size_t i = 0; i < x.size(); ++i)
        for (std::size_t j = 0; j < x[i].size(); ++j)
        {
            xPlus[i][j] += epsilon*direction[i][j];
            xMinus[i][j] -= epsilon*direction[i][j];
        }

    const auto residualPlus = residualAt(xPlus);
    const auto residualMinus = residualAt(xMinus);
    assembler.updateGridVariables(x);

    Scalar maxDeviation = 0.0, scale = 0.0;
    for (std::size_t i = 0; i < residual.size(); ++i)
        for (std::size_t j = 0; j < residual[i].size(); ++j)
        {
            using std::abs; using std::max;
            const auto finiteDifference = (residualPlus[i][j] - residualMinus[i][j])/(2.0*epsilon);
            maxDeviation = max(maxDeviation, abs(finiteDifference - jacobianTimesDirection[i][j]));
            scale = max(scale, abs(jacobianTimesDirection[i][j]));
        }

    return maxDeviation/scale;
}

} // end namespace Dumux::Test

#endif
