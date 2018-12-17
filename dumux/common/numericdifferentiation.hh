// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/****************************************************************************
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
 * \ingroup Common
 * \brief A class for numeric differentiation
 *
 */
#ifndef DUMUX_NUMERIC_DIFFERENTIATION_HH
#define DUMUX_NUMERIC_DIFFERENTIATION_HH

#include <cmath>
#include <limits>

namespace Dumux {

/*!
 * \ingroup Common
 * \brief A class for numeric differentiation with respect to a scalar parameter
 */
class NumericDifferentiation
{
public:

    /*!
     * \brief Computes the epsilon used for numeric differentiation
     * \param value The value of the variable with respect to which we are differentiating
     * \param baseEps The step width which we are using for differentiation
     */
    template<class Scalar>
    static Scalar epsilon(const Scalar value, const Scalar baseEps = 1e-10)
    {
        assert(std::numeric_limits<Scalar>::epsilon()*1e4 < baseEps);
        // the epsilon value used for the numeric differentiation is
        // now scaled by the absolute value of the primary variable...
        using std::abs;
        return baseEps*(abs(value) + 1.0);
    }

    /*!
     * \brief Computes the derivative of a function with repect to a function parameter
     * \note Overload using default epsilon computation
     */
    template<class Function, class Scalar, class FunctionEvalType>
    static void partialDerivative(const Function& function, Scalar x0,
                                  FunctionEvalType& derivative,
                                  const FunctionEvalType& fx0,
                                  const int numericDifferenceMethod = 1)
    { partialDerivative(function, x0, derivative, fx0, epsilon(x0), numericDifferenceMethod); }

    /*!
     * \brief Computes the derivative of a function with repect to a function parameter
     * \param function The function to derive
     * \param x0 The parameter at which the derivative is ought to be evaluated
     * \param derivative The partial derivative (output)
     * \param fx0 The result of the function evaluated at x0
     * \param eps The numeric epsilon used in the differentiation
     * \param numericDifferenceMethod The numeric difference method
     *        (1: foward differences (default), 0: central differences, -1: backward differences)
     */
    template<class Function, class Scalar, class FunctionEvalType>
    static void partialDerivative(const Function& function, Scalar x0,
                                  FunctionEvalType& derivative,
                                  const FunctionEvalType& fx0,
                                  const Scalar eps,
                                  const int numericDifferenceMethod = 1)
    {
        Scalar delta = 0.0;
        // we are using forward or central differences, i.e. we need to calculate f(x + \epsilon)
        if (numericDifferenceMethod >= 0)
        {
            delta += eps;
            // calculate the function evaluated with the deflected variable
            derivative = function(x0 + eps);
        }

        // we are using backward differences,
        // i.e. we don't need to calculate f(x + \epsilon)
        // we can recycle the (possibly cached) f(x)
        else derivative = fx0;

        // we are using backward or central differences,
        // i.e. we need to calculate f(x - \epsilon)
        if (numericDifferenceMethod <= 0)
        {
            delta += eps;
            // subtract the function evaluated with the deflected variable
            derivative -= function(x0 - eps);
        }

        // we are using forward differences,
        // i.e. we don't need to calculate f(x - \epsilon)
        // we can recycle the (possibly cached) f(x)
        else derivative -= fx0;

        // divide difference in residuals by the magnitude of the
        // deflections between the two function evaluation
        derivative /= delta;
    }
};

} // end namespace Dumux

#endif
