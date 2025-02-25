//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>

#include <iostream>
#include <cmath>
#include <random>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/math.hh>
#include <dumux/common/random.hh>

#include <dumux/nonlinear/leastsquares.hh>

namespace Dumux {

template<class F>
void testCurveFit(const Dune::BlockVector<double>& trueParameters, std::size_t numObservations, F&& function, const double noiseStddev = 1e-5, const double tol = 1e-3)
{
    // generate some observations
    const auto observationsX = Dumux::linspace(-1.0, 1.0, numObservations);

    // y = f(x) + noise
    // noise is drawn from a normal distribution with mean 0 and standard deviation noiseStddev
    std::mt19937 gen(0); // fixed seed for reproducibility
    SimpleNormalDistribution<double> noise(0.0, noiseStddev);

    auto observationsY = observationsX;
    std::transform(observationsX.begin(), observationsX.end(), observationsY.begin(),
                   [&](const double x) { return function(trueParameters, x) + noise(gen); });

    // residual function
    Dune::BlockVector<double> r(numObservations);
    const auto residual = [&](const auto& p)
    {
        for (std::size_t i = 0; i < numObservations; ++i)
            r[i] = observationsY[i] - function(p, observationsX[i]);
        return r;
    };

    // initial guess
    auto initialGuess(trueParameters);
    initialGuess = 1.0;

    auto optimizer = Optimization::makeNonlinearLeastSquaresSolver(residual, initialGuess, numObservations);
    auto params = initialGuess;
    const bool converged = optimizer->apply(params);
    if (!converged)
        DUNE_THROW(Dune::Exception, "Nonlinear least squares solver did not converge.");

    std::cout << "True parameters:\n" << trueParameters << std::endl;
    std::cout << "Estimated parameters:\n" << params << std::endl;

    // check if the estimated parameters are close to the true parameters
    for (std::size_t i = 0; i < 3; ++i)
        if (!Dune::FloatCmp::eq(params[i], trueParameters[i], tol))
            DUNE_THROW(Dune::Exception, "Estimated parameters are not close to the true parameters."
                << "Parameter: " << i << ", true: " << trueParameters[i] << ", estimated: " << params[i]);
}

} // end namespace Dumux

int main(int argc, char* argv[])
{
    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);

    // initialize parameters
    Dumux::Parameters::init(argc, argv);

    // simple curve fitting problems

    // f(x) = a*x + b * exp(c*x)
    {
        const auto trueParameters = Dune::BlockVector<double>({1.0, 2.0, 0.5});
        std::cout << "Test function f(x) = a*x + b * exp(c*x)" << std::endl;
        Dumux::testCurveFit(
            trueParameters, 50,
            [&](const auto& p, const double x)
            { return p[0]*x + p[1] * std::exp(p[2]*x); }
        );

        std::cout << "Test function f(x) = a*x + b * exp(c*x) with more noise" << std::endl;
        Dumux::testCurveFit(
            trueParameters, 50,
            [&](const auto& p, const double x)
            { return p[0]*x + p[1] * std::exp(p[2]*x); },
            0.01, 0.3
        );
    }

    // f(x) = a*x**3 + b*x**2 + c*x + d
    {
        std::cout << "Test function f(x) = a*x**3 + b*x**2 + c*x + d" << std::endl;
        const auto trueParameters = Dune::BlockVector<double>({1.0, -2.0, 3.0, 4.0});
        Dumux::testCurveFit(
            trueParameters, 50,
            [&](const auto& p, const double x)
            { return p[0]*std::pow(x, 3) + p[1]*std::pow(x, 2) + p[2]*x + p[3]; }
        );
    }

    // f(x) = a*sin(b*x + c)
    {
        std::cout << "Test function f(x) = a*sin(b*x + c)" << std::endl;
        const auto trueParameters = Dune::BlockVector<double>({1.0, 1.5, 0.7});
        Dumux::testCurveFit(
            trueParameters, 50,
            [&](const auto& p, const double x)
            { return p[0]*std::sin(p[1]*x + p[2]); }
        );
    }

    // quartic with large coeffs
    {
        std::cout << "Test function f(x) = a*x**4 + b*x**3 + c*x**2 + d*x + e" << std::endl;
        const auto trueParameters = Dune::BlockVector<double>({1e5, -2e4, 3e3, 4e2, 5e1});
        Dumux::testCurveFit(
            trueParameters, 50,
            [&](const auto& p, const double x)
            { return p[0]*std::pow(x, 4) + p[1]*std::pow(x, 3) + p[2]*std::pow(x, 2) + p[3]*x + p[4]; }
        );
    }

    return 0;
}
