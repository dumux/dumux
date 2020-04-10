// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Nonlinear
 * \brief Root finding algorithms for scalar functions
 */
#ifndef DUMUX_COMMON_SCALAR_ROOT_FINDING_HH
#define DUMUX_COMMON_SCALAR_ROOT_FINDING_HH

#include <cmath>
#include <limits>
#include <type_traits>

#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numericdifferentiation.hh>

namespace Dumux {

/*!
 * \ingroup Nonlinear
 * \brief Newton's root finding algorithm for scalar functions (secant method)
 * \param xOld initial guess
 * \param residual Residual function
 * \param derivative Derivative of the residual
 * \param tol Relative shift tolerance
 * \param maxIter Maximum number of iterations
 */
template<class Scalar, class ResFunc, class DerivFunc,
         typename std::enable_if_t<std::is_invocable_r_v<Scalar, ResFunc, Scalar> &&
                                   std::is_invocable_r_v<Scalar, DerivFunc, Scalar>>...>
Scalar findScalarRootNewton(Scalar xOld, const ResFunc& residual, const DerivFunc& derivative,
                            const Scalar tol = 1e-13, const int maxIter = 200)
{
    Scalar xNew = xOld;
    Scalar r = residual(xNew);

    int n = maxIter;
    Scalar relativeShift = std::numeric_limits<Scalar>::max();
    while (relativeShift > tol && n > 0)
    {
        xNew = xOld - r/derivative(xOld);
        r = residual(xNew);

        using std::abs; using std::max;
        relativeShift = abs(xOld-xNew)/max(abs(xOld), abs(xNew));
        xOld = xNew;
        n--;
    }

    using std::isfinite;
    if (!isfinite(r))
        DUNE_THROW(NumericalProblem, "Residual is not finite: " << r << " after " << maxIter - n << " iterations!");

    if (relativeShift > tol)
        DUNE_THROW(NumericalProblem, "Scalar newton solver did not converge after " << maxIter << " iterations!");

    return xNew;
}

/*!
 * \ingroup Nonlinear
 * \brief Newton's root finding algorithm for scalar functions (secant method)
 * \note The derivative is numerically computed. If the derivative is know use signature with derivative function.
 * \param xOld initial guess
 * \param residual Residual function
 * \param tol Relative shift tolerance
 * \param maxIter Maximum number of iterations
 */
template<class Scalar, class ResFunc,
          typename std::enable_if_t<std::is_invocable_r_v<Scalar, ResFunc, Scalar>>...>
Scalar findScalarRootNewton(Scalar xOld, const ResFunc& residual,
                            const Scalar tol = 1e-13, const int maxIter = 200)
{
    const auto eps = NumericDifferentiation::epsilon(xOld);
    auto derivative = [&](const auto x){ return (residual(x + eps)-residual(x))/eps; };
    return findScalarRootNewton(xOld, residual, derivative, tol, maxIter);
}

/*!
 * \ingroup Nonlinear
 * \brief Brent's root finding algorithm for scalar functions
 * \note Modified from pseudo-code on wikipedia: https://en.wikipedia.org/wiki/Brent%27s_method
 * \note See also R.P. Brent "An algorithm with guaranteed convergence for finding a zero of a function", The Computer Journal (1971).
 * \note This is usually more robust than Newton's method
 * \param a Lower bound
 * \param b Upper bound
 * \param residual Residual function
 * \param tol Relative shift tolerance
 * \param maxIter Maximum number of iterations
 */
template<class Scalar, class ResFunc,
         typename std::enable_if_t<std::is_invocable_r_v<Scalar, ResFunc, Scalar>>...>
Scalar findScalarRootBrent(Scalar a, Scalar b, const ResFunc& residual,
                           const Scalar tol = 1e-13, const int maxIter = 200)
{
    // precompute the residuals (minimize function evaluations)
    Scalar fa = residual(a);
    Scalar fb = residual(b);
    Scalar fs = 0.0;

    // check if the root is inside the given interval
    using std::signbit;
    if (!signbit(fa*fb))
        DUNE_THROW(NumericalProblem, "Brent's algorithm failed: [a,b] does not contain any, or no uniquely defined root!");

    // sort bounds
    using std::abs; using std::swap;
    if (abs(fa) < abs(fb))
    {
        swap(a, b);
        swap(fa, fb);
    }

    Scalar c = a;
    Scalar fc = fa;
    Scalar d = 0.0;
    Scalar s = 0.0;
    bool flag = true;

    for (int i = 0; i < maxIter; ++i)
    {
        // stopping criterion
        using std::max;
        if (abs(b-a) < tol*max(abs(a), abs(b)))
            return b;

        // inverse quadratic interpolation
        if (fa != fc && fb != fc)
        {
            const auto fab = fa-fb;
            const auto fac = fa-fc;
            const auto fbc = fb-fc;
            s = a*fb*fc/(fab*fac) - b*fa*fc/(fab*fbc) + c*fa*fb/(fac*fbc);
        }

        // secant method
        else
        {
            s = b - fb*(b-a)/(fb-fa);
        }

        // bisection method
        if ( (s < (3*a + b)*0.25 || s > b)
             || (flag && abs(s-b) >= abs(b-c)*0.5)
             || (!flag && abs(s-b) >= abs(c-d)*0.5)
             || (flag && abs(b-c) < tol*max(abs(b), abs(c)))
             || (!flag && abs(c-d) < tol*max(abs(c), abs(d))) )
        {
            s = (a+b)*0.5;
            flag = true;
        }
        else
            flag = false;

        // compute residual at new guess
        fs = residual(s);
        d = c;
        c = b;
        fc = fb;

        // check on which side of the root s falls
        if (fa*fs < 0.0)
        {
            b = s;
            fb = fs;
        }
        else
        {
            a = s;
            fa = fs;
        }

        // sort if necessary
        if (abs(fa) < abs(fb))
        {
            swap(a, b);
            swap(fa, fb);
        }
    }

    DUNE_THROW(NumericalProblem, "Brent's algorithm didn't converge after " << maxIter << " iterations!");
}

} // end namespace Dumux

#endif
