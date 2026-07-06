// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PhaseFieldTests
 * \brief Test for the phase-field free-energy potentials.
 */
#include <config.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string_view>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dumux/common/numericdifferentiation.hh>

#include <dumux/phasefield/freeenergy/freeenergy.hh>

namespace {

template<class F, class D>
void checkDerivative(std::string_view name, const std::vector<double>& values, const F& f, const D& df)
{
    // absolute tolerance: these potentials are O(1) in magnitude and some have
    // exact zero-crossings (e.g. at c=1/2), where a relative comparison against
    // the O(1e-9) numerical-differentiation noise floor would spuriously fail
    static constexpr double eps = 1e-5;
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

/*!
 * \brief Checks a convex-concave (Eyre) splitting of a potential: that the
 *        two parts recombine exactly into the full derivative, and that both
 *        parts are genuinely convex (non-decreasing derivative).
 */
template<class F>
void checkConvexConcaveSplit(std::string_view name, const std::vector<double>& values, const F& f)
{
    static constexpr double eps = 1e-10;
    static constexpr double numEps = 1e-8;

    for (const auto c : values)
    {
        const auto split = f.convexDerivative(c) - f.concaveDerivative(c);
        using std::abs;
        if (abs(split - f.derivative(c)) > eps)
            DUNE_THROW(Dune::Exception, "Convex-concave split does not recombine into derivative() for "
                       << name << ": " << split << " != " << f.derivative(c) << " at c=" << c);
    }

    // convexity <=> non-decreasing derivative; check via numeric differentiation
    // (a small negative tolerance accommodates numerical-differentiation noise)
    static constexpr double convexityTol = -1e-6;
    for (const auto c : values)
    {
        double convexSecondDerivative = 0.0;
        Dumux::NumericDifferentiation::partialDerivative(
            [&](auto c){ return f.convexDerivative(c); }, c, convexSecondDerivative, f.convexDerivative(c), numEps, 0);
        if (convexSecondDerivative < convexityTol)
            DUNE_THROW(Dune::Exception, "Convex part of " << name << " is not convex at c=" << c
                       << ": f_vex''=" << convexSecondDerivative);

        double concaveSecondDerivative = 0.0;
        Dumux::NumericDifferentiation::partialDerivative(
            [&](auto c){ return f.concaveDerivative(c); }, c, concaveSecondDerivative, f.concaveDerivative(c), numEps, 0);
        if (concaveSecondDerivative < convexityTol)
            DUNE_THROW(Dune::Exception, "Concave part of " << name << " is not convex at c=" << c
                       << ": f_cave''=" << concaveSecondDerivative);
    }
}

} // end anonymous namespace

int main()
{
    using namespace Dumux::PhaseField::FreeEnergy;
    using std::abs;

    const auto c = linspace(0.02, 0.98, 49);

    {
        DoubleWell<double> f(2.5);
        checkDerivative("DoubleWell::derivative", c, [&](auto c){ return f.value(c); }, [&](auto c){ return f.derivative(c); });
        checkDerivative("DoubleWell::secondDerivative", c, [&](auto c){ return f.derivative(c); }, [&](auto c){ return f.secondDerivative(c); });

        if (abs(f.value(0.0)) > 1e-14 || abs(f.value(1.0)) > 1e-14)
            DUNE_THROW(Dune::Exception, "DoubleWell should vanish at c=0 and c=1");
        if (abs(f.derivative(0.0)) > 1e-14 || abs(f.derivative(1.0)) > 1e-14)
            DUNE_THROW(Dune::Exception, "DoubleWell derivative should vanish at c=0 and c=1 (minima)");
        if (abs(f.derivative(0.5)) > 1e-14)
            DUNE_THROW(Dune::Exception, "DoubleWell derivative should vanish at c=1/2 (local maximum)");

        checkConvexConcaveSplit("DoubleWell", c, f);
    }

    {
        DoubleObstacle<double> f(1.5);
        checkDerivative("DoubleObstacle::derivative", c, [&](auto c){ return f.value(c); }, [&](auto c){ return f.derivative(c); });
        checkDerivative("DoubleObstacle::secondDerivative", c, [&](auto c){ return f.derivative(c); }, [&](auto c){ return f.secondDerivative(c); });

        if (abs(f.derivative(0.5)) > 1e-14)
            DUNE_THROW(Dune::Exception, "DoubleObstacle derivative should vanish at c=1/2");
    }

    {
        Logarithmic<double> f(0.6, 1.0);
        checkDerivative("Logarithmic::derivative", c, [&](auto c){ return f.value(c); }, [&](auto c){ return f.derivative(c); });
        checkDerivative("Logarithmic::secondDerivative", c, [&](auto c){ return f.derivative(c); }, [&](auto c){ return f.secondDerivative(c); });

        if (abs(f.derivative(0.5)) > 1e-14)
            DUNE_THROW(Dune::Exception, "Logarithmic derivative should vanish at c=1/2 (symmetric potential)");
    }

    {
        // sample on [-1,1] instead of [0,1] -- this potential uses the symmetric range
        const auto phi = linspace(-0.98, 0.98, 49);
        SymmetricDoubleWell<double> f(2.5);
        checkDerivative("SymmetricDoubleWell::derivative", phi, [&](auto phi){ return f.value(phi); }, [&](auto phi){ return f.derivative(phi); });
        checkDerivative("SymmetricDoubleWell::secondDerivative", phi, [&](auto phi){ return f.derivative(phi); }, [&](auto phi){ return f.secondDerivative(phi); });

        if (abs(f.value(-1.0)) > 1e-14 || abs(f.value(1.0)) > 1e-14)
            DUNE_THROW(Dune::Exception, "SymmetricDoubleWell should vanish at phi=-1 and phi=1");
        if (abs(f.derivative(-1.0)) > 1e-14 || abs(f.derivative(1.0)) > 1e-14)
            DUNE_THROW(Dune::Exception, "SymmetricDoubleWell derivative should vanish at phi=-1 and phi=1 (minima)");
        if (abs(f.derivative(0.0)) > 1e-14)
            DUNE_THROW(Dune::Exception, "SymmetricDoubleWell derivative should vanish at phi=0 (local maximum)");

        checkConvexConcaveSplit("SymmetricDoubleWell", phi, f);
    }

    std::cout << "All free-energy potential tests passed." << std::endl;
    return 0;
}
