// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PhaseFieldTests
 * \brief Test for the phase-field mobility functions.
 */
#include <config.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string_view>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dumux/common/numericdifferentiation.hh>

#include <dumux/phasefield/mobility/mobility.hh>

namespace {

template<class F, class D>
void checkDerivative(std::string_view name, const std::vector<double>& values, const F& f, const D& df)
{
    // absolute tolerance: these mobilities are O(1e-4) in magnitude and
    // Degenerate has an exact zero-crossing at c=1/2, where a relative
    // comparison against the numerical-differentiation noise floor would
    // spuriously fail
    static constexpr double eps = 1e-9;
    static constexpr double numEps = 1e-8;
    for (const auto c : values)
    {
        const double analytic = df(c);
        double numeric = 0.0;
        Dumux::NumericDifferentiation::partialDerivative(f, c, numeric, f(c), numEps, 0 /*central differences*/);

        using std::abs;
        if (abs(analytic - numeric) > eps)
            DUNE_THROW(Dune::Exception, "Derivative mismatch for " << name << ": "
                       << std::setprecision(10) << analytic << " != " << numeric << " at c=" << c);
    }
}

std::vector<double> linspace(double a, double b, std::size_t n)
{
    std::vector<double> values(n);
    for (std::size_t i = 0; i < n; ++i)
        values[i] = a + (b-a)*double(i)/double(n-1);
    return values;
}

} // end anonymous namespace

int main()
{
    using namespace Dumux::PhaseField::Mobility;
    using std::abs;

    const auto c = linspace(0.02, 0.98, 49);

    {
        Constant<double> m(1e-4);
        checkDerivative("Constant::derivative", c, [&](auto c){ return m.value(c); }, [&](auto c){ return m.derivative(c); });

        for (const auto ci : c)
            if (abs(m.value(ci) - 1e-4) > 1e-14)
                DUNE_THROW(Dune::Exception, "Constant mobility must equal its value everywhere");
    }

    {
        Degenerate<double> m(1e-4);
        checkDerivative("Degenerate::derivative", c, [&](auto c){ return m.value(c); }, [&](auto c){ return m.derivative(c); });

        if (abs(m.value(0.0)) > 1e-14 || abs(m.value(1.0)) > 1e-14)
            DUNE_THROW(Dune::Exception, "Degenerate mobility should vanish at c=0 and c=1");

        for (const auto ci : c)
            if (m.value(ci) < 0.0)
                DUNE_THROW(Dune::Exception, "Degenerate mobility should be non-negative on [0,1]");
    }

    std::cout << "All mobility tests passed." << std::endl;
    return 0;
}
