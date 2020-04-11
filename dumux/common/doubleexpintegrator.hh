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
 *                                                                           *
 *   This version is modified after the original version by John D. Cook     *
 *   see https://www.codeproject.com/                                        *
 *       Articles/31550/Fast-Numerical-Integration                           *
 *   which is licensed under BSD-2-clause, which reads as follows:           *
 *   Copyright John D. Cook                                                  *
 *                                                                           *
 *   Redistribution and use in source and binary forms, with or without      *
 *   modification, are permitted provided that the following                 *
 *   conditions are met:                                                     *
 *   1. Redistributions of source code must retain the above                 *
 *      copyright notice, this list of conditions                            *
        and the following disclaimer.                                        *
 *   2. Redistributions in binary form must reproduce the above              *
 *      copyright notice, this list of conditions                            *
 *      and the following disclaimer                                         *
 *   in the documentation and/or other materials                             *
 *   provided with the distribution.                                         *
 *                                                                           *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS     *
 *   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT       *
 *   LIMITED TO, THE IMPLIED  WARRANTIES OF MERCHANTABILITY AND FITNESS      *
 *   FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE          *
 *   COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,              *
 *   INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL              *
 *   DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE       *
 *   GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS           *
 *   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,            *
 *   WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE    *
 *   OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,       *
 *   EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                      *
 *                                                                           *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Common
 * \brief A double exponential integrator
 */
#ifndef DUMUX_COMMON_DOUBLEEXP_INTEGRATOR_HH
#define DUMUX_COMMON_DOUBLEEXP_INTEGRATOR_HH

#include <cmath>
#include <limits>
#include <type_traits>

#include <dumux/common/doubleexpintegrationconstants.hh>

namespace Dumux {

/*!
 * \brief Numerical integration in one dimension using the double exponential method of M. Mori.
 */
template<class Scalar>
class DoubleExponentialIntegrator
{
public:
    /*!
     * \brief Integrate an analytic function over a finite interval
     * \param f the integrand (invocable with a single scalar)
     * \param a lower limit of integration
     * \param b upper limit of integration
     * \param targetAbsoluteError desired bound on error
     * \param numFunctionEvaluations number of function evaluations used
     * \param errorEstimate estimate for error in integration
     * \return The value of the integral
     */
    template<class Function,
             typename std::enable_if_t<std::is_invocable_r_v<Scalar, Function, Scalar>>...>
    static Scalar integrate(const Function& f,
                            const Scalar a,
                            const Scalar b,
                            const Scalar targetAbsoluteError,
                            int& numFunctionEvaluations,
                            Scalar& errorEstimate)
    {
        // Apply the linear change of variables x = ct + d
        // $$\int_a^b f(x) dx = c \int_{-1}^1 f( ct + d ) dt$$
        // c = (b-a)/2, d = (a+b)/2

        const Scalar c = 0.5*(b - a);
        const Scalar d = 0.5*(a + b);
        return integrateCore_(f, c, d, targetAbsoluteError, numFunctionEvaluations, errorEstimate,
                              doubleExponentialIntegrationAbcissas, doubleExponentialIntegrationWeights);
    }

    /*!
     * \brief Integrate an analytic function over a finite interval.
     * \note This version overloaded to not require arguments passed in for
     *       function evaluation counts or error estimates.
     * \param f the integrand (invocable with a single scalar)
     * \param a lower integral bound
     * \param b upper integral bound
     * \param targetAbsoluteError desired absolute error in the result
     * \return The value of the integral.
     */
    template<class Function,
             typename std::enable_if_t<std::is_invocable_r_v<Scalar, Function, Scalar>>...>
    static Scalar integrate(const Function& f,
                            const Scalar a,
                            const Scalar b,
                            const Scalar targetAbsoluteError)
    {
        int numFunctionEvaluations;
        Scalar errorEstimate;
        return integrate(f, a, b, targetAbsoluteError, numFunctionEvaluations, errorEstimate);
    }

private:
    // Integrate f(cx + d) with the given integration constants
    template<class Function>
    static Scalar integrateCore_(const Function& f,
                                 const Scalar c,   // slope of change of variables
                                 const Scalar d,   // intercept of change of variables
                                 Scalar targetAbsoluteError,
                                 int& numFunctionEvaluations,
                                 Scalar& errorEstimate,
                                 const double* abcissas,
                                 const double* weights)
    {
        targetAbsoluteError /= c;

        // Offsets to where each level's integration constants start.
        // The last element is not a beginning but an end.
        static const int offsets[] = {1, 4, 7, 13, 25, 49, 97, 193};
        static const int numLevels = sizeof(offsets)/sizeof(int) - 1;

        Scalar newContribution = 0.0;
        Scalar integral = 0.0;
        Scalar h = 1.0;
        errorEstimate = std::numeric_limits<Scalar>::max();
        Scalar previousDelta, currentDelta = std::numeric_limits<Scalar>::max();

        integral = f(c*abcissas[0] + d) * weights[0];
        int i;
        for (i = offsets[0]; i != offsets[1]; ++i)
            integral += weights[i]*(f(c*abcissas[i] + d) + f(-c*abcissas[i] + d));

        for (int level = 1; level != numLevels; ++level)
        {
            h *= 0.5;
            newContribution = 0.0;
            for (i = offsets[level]; i != offsets[level+1]; ++i)
                newContribution += weights[i]*(f(c*abcissas[i] + d) + f(-c*abcissas[i] + d));
            newContribution *= h;

            // difference in consecutive integral estimates
            previousDelta = currentDelta;
            using std::abs;
            currentDelta = abs(0.5*integral - newContribution);
            integral = 0.5*integral + newContribution;

            // Once convergence kicks in, error is approximately squared at each step.
            // Determine whether we're in the convergent region by looking at the trend in the error.
            if (level == 1)
                continue; // previousDelta meaningless, so cannot check convergence.

            // Exact comparison with zero is harmless here.  Could possibly be replaced with
            // a small positive upper limit on the size of currentDelta, but determining
            // that upper limit would be difficult.  At worse, the loop is executed more
            // times than necessary.  But no infinite loop can result since there is
            // an upper bound on the loop variable.
            if (currentDelta == 0.0)
                break;

            using std::log;
            const Scalar rate = log( currentDelta )/log( previousDelta );  // previousDelta != 0 or would have been kicked out previously

            if (rate > 1.9 && rate < 2.1)
            {
                // If convergence theory applied perfectly, r would be 2 in the convergence region.
                // r close to 2 is good enough. We expect the difference between this integral estimate
                // and the next one to be roughly delta^2.
                errorEstimate = currentDelta*currentDelta;
            }
            else
            {
                // Not in the convergence region.  Assume only that error is decreasing.
                errorEstimate = currentDelta;
            }

            if (errorEstimate < 0.1*targetAbsoluteError)
                break;
        }

        numFunctionEvaluations = 2*i - 1;
        errorEstimate *= c;
        return c*integral;
    }
};

} // end namespace Dumux

#endif
