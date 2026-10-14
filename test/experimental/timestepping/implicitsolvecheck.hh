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
 * For a given solution, the residual assembled for the weights (0, wM, wR) and a
 * constant c has to equal wM M + wR R + c in all unconstrained equations, where
 * M and R are the operators returned by assembleOperators, and has to be independent
 * of c in the constrained equations.
 * \return the number of constrained equations
 */
template<class Assembler, class SolutionVector, class Scalar>
std::size_t checkImplicitSolveContract(Assembler& assembler, const SolutionVector& x, const Scalar time,
                                       const Scalar storageWeight, const Scalar spatialWeight)
{
    using Residual = typename Assembler::ResidualType;
    using SolveParams = Experimental::ImplicitSolveParams<Scalar, Residual>;
    const Experimental::TimeLevel<Scalar> timeLevel(time);

    Residual storage, spatial;
    assembler.assembleOperators(x, timeLevel, {.storage = &storage, .spatial = &spatial});

    auto constant = std::make_shared<Residual>(storage);
    std::size_t counter = 0;
    Detail::forEachEntry(*constant, storage, storage, storage, [&](auto& c, auto&, auto&, auto&)
    { c = std::sin(1.0 + 0.37*(counter++)); });

    assembler.setImplicitSolve(x, SolveParams{timeLevel, {0.0, storageWeight, spatialWeight}, nullptr});
    assembler.assembleJacobianAndResidual(x);
    Residual withoutConstant = assembler.residual();

    assembler.setImplicitSolve(x, SolveParams{timeLevel, {0.0, storageWeight, spatialWeight}, constant});
    assembler.assembleJacobianAndResidual(x);
    Residual withConstant = assembler.residual();

    Scalar scale = 0.0;
    Detail::forEachEntry(withoutConstant, storage, spatial, withConstant, [&](auto& r, auto& m, auto& s, auto&)
    {
        using std::abs; using std::max;
        scale = max({scale, abs(r), abs(storageWeight*m), abs(spatialWeight*s)});
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

} // end namespace Dumux::Test

#endif
