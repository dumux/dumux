// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Experimental
 * \brief Roots of polynomials with complex coefficients for the time stepping tests
 */
#ifndef DUMUX_TEST_TIMESTEPPING_POLYNOMIAL_ROOTS_HH
#define DUMUX_TEST_TIMESTEPPING_POLYNOMIAL_ROOTS_HH

#include <cmath>
#include <vector>
#include <complex>
#include <algorithm>

namespace Dumux::Test {

//! The roots of sum_j c_j zeta^{n-j} (Durand-Kerner iteration)
inline std::vector<std::complex<double>> roots(const std::vector<std::complex<double>>& coefficients)
{
    const auto n = coefficients.size() - 1;
    std::vector<std::complex<double>> z(n);
    for (std::size_t i = 0; i < n; ++i)
        z[i] = std::pow(std::complex<double>{0.4, 0.9}, double(i));

    const auto evaluate = [&](const std::complex<double> x)
    {
        std::complex<double> value = 0.0;
        for (const auto& c : coefficients)
            value = value*x + c/coefficients[0];
        return value;
    };

    for (int iteration = 0; iteration < 1000; ++iteration)
    {
        double change = 0.0;
        for (std::size_t i = 0; i < n; ++i)
        {
            std::complex<double> denominator = 1.0;
            for (std::size_t j = 0; j < n; ++j)
                if (j != i)
                    denominator *= z[i] - z[j];
            const auto delta = evaluate(z[i])/denominator;
            z[i] -= delta;
            change = std::max(change, std::abs(delta));
        }
        if (change < 1e-15)
            break;
    }

    return z;
}

} // end namespace Dumux::Test

#endif
