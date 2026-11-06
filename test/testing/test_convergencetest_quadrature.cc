// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Errors of the midpoint rule (second order) and the left Riemann sum (first order)
 *        for the integral of exp(x) on [0, 1], and of the midpoint rule for the integral of x,
 *        which it integrates exactly, written for the convergence test driver
 *
 * Usage: test_convergencetest_quadrature <name> <number of intervals> [nan]
 */
#include <config.h>

#include <cmath>
#include <limits>
#include <string>

#include <test/testing/convergencetest.hh>

int main(int argc, char** argv)
{
    if (argc < 3)
        return 1;

    const std::string name = argv[1];
    const int numIntervals = std::stoi(argv[2]);
    if (numIntervals < 1)
        return 2;

    const double h = 1.0/numIntervals;
    double midpoint = 0.0;
    double leftSum = 0.0;
    double midpointLinear = 0.0;
    for (int i = 0; i < numIntervals; ++i)
    {
        midpoint += h*std::exp((i + 0.5)*h);
        leftSum += h*std::exp(i*h);
        midpointLinear += h*(i + 0.5)*h;
    }

    const double exact = std::exp(1.0) - 1.0;
    const bool writeNaN = argc > 3 && std::string(argv[3]) == "nan";
    Dumux::Test::writeConvergenceTestErrors(name, h, {
        {"midpoint", writeNaN ? std::numeric_limits<double>::quiet_NaN() : std::abs(midpoint - exact)},
        {"leftsum", std::abs(leftSum - exact)},
        {"linear", std::abs(midpointLinear - 0.5)}
    });

    return 0;
}
