//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>

#include <dumux/common/numericdifferentiation.hh>

#include <cmath>
#include <iostream>
#include <limits>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>

#include <dumux/io/format.hh>

int main(int argc, char* argv[])
{
    using namespace Dumux;

    // There is two sources of error in the numerical differentiation:
    // 1. The truncation error, which is the error introduced by the approximation of the derivative
    // 2. The round-off and representation errors, which are the errors introduced by the finite precision of the computer
    // There is a trade-off between the two errors, as the truncation error decreases with a smaller step size
    // while the round-off error increases with a smaller step size.

    // For the higher-order methods, larger step sizes can/should be used, as the truncation error is smaller
    // and this then avoids the round-off error to dominate the total error.
    // For a central difference scheme, the step size eps = √ε where ε: machine precision epsilon, is optimal.

    const auto testNumDiff = [](const auto& func, const auto& dfunc,
                                const double x0, const double eps, const int numDiffMethod,
                                const double tol)
    {
        // scalar type //////
        const double fx0 = func(x0);
        double dfx0 = 0.0;
        NumericDifferentiation::partialDerivative(func, x0, dfx0, fx0, eps, numDiffMethod);
        const double dfx0Exact = dfunc(x0);

        std::cout << Fmt::format("Scalar test: f'(1) = {:.15e}; exact: {:.15e}\n", dfx0, dfx0Exact)
                  << Fmt::format("-- eps = {:.15e}, tol = {:.15e}, method = {}\n", eps, tol, numDiffMethod);

        if (!Dune::FloatCmp::eq(dfx0, dfx0Exact, tol))
        {
            std::cout << "--> FAILED\n";
            return 1;
        }

        // vector type //////
        const auto funcVec = [&](double x) -> Dune::FieldVector<double, 2> { return { func(x), func(x) }; };
        const auto dfuncVec = [&](double x) -> Dune::FieldVector<double, 2> { return { dfunc(x), dfunc(x) }; };
        const Dune::FieldVector<double, 2> fx0Vec = funcVec(x0);

        Dune::FieldVector<double, 2> dfx0Vec = Dune::FieldVector<double, 2>(0.0);
        NumericDifferentiation::partialDerivative(funcVec, x0, dfx0Vec, fx0Vec, eps, numDiffMethod);
        const Dune::FieldVector<double, 2> dfx0ExactVec = dfuncVec(x0);

        std::cout << Fmt::format("Vector test: f'(1) = {:.15e}, {:.15e}; exact: {:.15e}, {:.15e}\n", dfx0Vec[0], dfx0Vec[1], dfx0ExactVec[0], dfx0ExactVec[1])
                  << Fmt::format("-- eps = {:.15e}, tol = {:.15e}, method = {}\n", eps, tol, numDiffMethod);

        if (!Dune::FloatCmp::eq(dfx0Vec, dfx0ExactVec, tol))
        {
            std::cout << "--> FAILED\n";
            return 1;
        }

        return 0;
    };

    int failedTest = 0;

    {
        const auto func = [](double x){ return 1/24.0*x*x*x*x + 1/6.0*x*x*x + 1/2.0*x*x + x; };
        const auto dfunc = [](double x){ return 1/6.0*x*x*x + 1/2.0*x*x + x + 1; };
        const double x0 = 1.0;

        const double baseEps = std::sqrt(std::numeric_limits<double>::epsilon());
        const auto tests = std::vector<std::tuple<int, double, double>>({
            { 1, baseEps, 100*baseEps },
            { -1, baseEps, 100*baseEps },
            { 0, baseEps, baseEps }, // optimal step size for central difference is sqrt(eps)
            { 5, 1e6*baseEps, baseEps*baseEps } // for large enough step size the error is zero as 5th derivative of func is zero
        });

        for (const auto& [numDiffMethod, eps, tol] : tests)
            failedTest += testNumDiff(func, dfunc, x0, eps, numDiffMethod, tol);
    }

    if (failedTest > 0)
    {
        std::cout << Fmt::format("{} failed tests\n", failedTest);
        return 1;
    }

    std::cout << "All tests passed" << std::endl;

    return 0;
}
