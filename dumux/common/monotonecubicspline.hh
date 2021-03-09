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
 * \ingroup Common
 * \brief A monotone cubic spline
 */
#ifndef DUMUX_COMMON_MONOTONE_CUBIC_SPLINE_HH
#define DUMUX_COMMON_MONOTONE_CUBIC_SPLINE_HH

#include <cmath>
#include <cassert>
#include <iterator>
#include <vector>
#include <algorithm>

#include <dune/common/float_cmp.hh>

// Hermite basis functions
#include <dumux/common/cubicsplinehermitebasis.hh>

// for inversion
#include <dumux/nonlinear/findscalarroot.hh>

namespace Dumux {

/*!
 * \ingroup Common
 * \brief A monotone cubic spline
 * \note Construction after Fritsch & Butland (1984) (see https://doi.org/10.1137/0905021)
 * \note The resulting interpolation is globally monotone but only C^1
 */
template<class Scalar = double>
class MonotoneCubicSpline
{
    using Basis = CubicSplineHermiteBasis<Scalar>;
public:
    /*!
     * \brief Default constructor
     */
    MonotoneCubicSpline() = default;

    /*!
     * \brief Construct a monotone cubic spline from the control points (x[i], y[i])
     * \note if the data set is monotone, monotonicity is preserved
     * \param x a vector of x-coordinates
     * \param y a vector of y-coordinates
     */
    MonotoneCubicSpline(const std::vector<Scalar>& x, const std::vector<Scalar>& y)
    {
        updatePoints(x, y);
    }

    /*!
     * \brief Create a monotone cubic spline from the control points (x[i], y[i])
     * \param x a vector of x-coordinates
     * \param y a vector of y-coordinates
     */
    void updatePoints(const std::vector<Scalar>& x, const std::vector<Scalar>& y)
    {
        // check some requirements
        assert (x.size() == y.size());
        assert (x.size() >=2);
        assert (std::is_sorted(x.begin(), x.end()) || std::is_sorted(x.rbegin(), x.rend()));

        // save a copy of the control points
        x_ = x;
        y_ = y;

        // the number of control points
        numPoints_ = x.size();

        // whether we are increasing
        increasingX_ = x_.back() > x_.front();
        increasingY_ = y_.back() > y_.front();

        // the slope at every control point
        m_.resize(numPoints_);

        // compute slopes (see Fritsch & Butland (1984), Eq. (5))
        Scalar deltaX = (x[1]-x[0]);
        Scalar secant = m_.front() = (y[1]-y[0])/deltaX;
        Scalar prevDeltaX = deltaX;
        Scalar prevSecant = secant;
        for (int i = 1; i < numPoints_-1; ++i, prevSecant = secant, prevDeltaX = deltaX)
        {
            deltaX = (x[i+1]-x[i]);
            secant = (y[i+1]-y[i])/deltaX;
            const auto alpha = (prevDeltaX + 2*deltaX)/(3*(prevDeltaX + deltaX));
            m_[i] = prevSecant*secant > 0.0 ? prevSecant*secant/(alpha*secant + (1.0-alpha)*prevSecant) : 0.0;
        }
        m_.back() = secant;
    }

    /*!
     * \brief Evaluate the y value at a given x value
     * \param x the x-coordinate
     * \note We extrapolate linearly if out of bounds
     */
    Scalar eval(const Scalar x) const
    {
        if ((x <= x_.front() && increasingX_) || (x >= x_.front() && !increasingX_))
            return y_.front() + m_.front()*(x - x_.front());
        else if ((x > x_.back() && increasingX_) || (x < x_.back() && !increasingX_))
            return y_.back() + m_.back()*(x - x_.back());

        return eval_(x);
    }

    /*!
     * \brief Evaluate the first derivative dy/dx at a given x value
     * \param x the x-coordinate
     * \note We extrapolate linearly if out of bounds
     */
    Scalar evalDerivative(const Scalar x) const
    {
        if ((x <= x_.front() && increasingX_) || (x >= x_.front() && !increasingX_))
            return m_.front();
        else if ((x > x_.back() && increasingX_) || (x < x_.back() && !increasingX_))
            return m_.back();

        return evalDerivative_(x);
    }

    /*!
     * \brief Evaluate the inverse function
     * \param y the y-coordinate
     * \note We extrapolate linearly if out of bounds
     * \note Throws exception if inverse could not be found (e.g. not unique)
     */
    Scalar evalInverse(const Scalar y) const
    {
        if ((y <= y_.front() && increasingY_) || (y >= y_.front() && !increasingY_))
            return x_.front() + (y - y_.front())/m_.front();
        else if ((y > y_.back() && increasingY_) || (y < y_.back() && !increasingY_))
            return x_.back() + (y - y_.back())/m_.back();

        return evalInverse_(y);
    }

private:
    Scalar eval_(const Scalar x) const
    {
        // interpolate parametrization parameter t in [0,1]
        const auto lookUpIndex = lookUpIndex_(x_, x, increasingX_);
        const auto h = (x_[lookUpIndex] - x_[lookUpIndex-1]);
        const auto t = (x - x_[lookUpIndex-1])/h;
        return y_[lookUpIndex-1]*Basis::h00(t) + h*m_[lookUpIndex-1]*Basis::h10(t)
               + y_[lookUpIndex]*Basis::h01(t) + h*m_[lookUpIndex]*Basis::h11(t);
    }

    Scalar evalDerivative_(const Scalar x) const
    {
        // interpolate parametrization parameter t in [0,1]
        const auto lookUpIndex = lookUpIndex_(x_, x, increasingX_);
        const auto h = (x_[lookUpIndex] - x_[lookUpIndex-1]);
        const auto t = (x - x_[lookUpIndex-1])/h;
        const auto dtdx = 1.0/h;
        return y_[lookUpIndex-1]*Basis::dh00(t)*dtdx + m_[lookUpIndex-1]*Basis::dh10(t)
               + y_[lookUpIndex]*Basis::dh01(t)*dtdx + m_[lookUpIndex]*Basis::dh11(t);
    }

    Scalar evalInverse_(const Scalar y) const
    {
        const auto lookUpIndex = lookUpIndex_(y_, y, increasingY_);
        auto localPolynomial = [&](const auto x) {
            // interpolate parametrization parameter t in [0,1]
            const auto h = (x_[lookUpIndex] - x_[lookUpIndex-1]);
            const auto t = (x - x_[lookUpIndex-1])/h;
            return y - (y_[lookUpIndex-1]*Basis::h00(t) + h*m_[lookUpIndex-1]*Basis::h10(t)
                   + y_[lookUpIndex]*Basis::h01(t) + h*m_[lookUpIndex]*Basis::h11(t));
        };

        // use an epsilon for the bracket
        const auto eps = (x_[lookUpIndex]-x_[lookUpIndex-1])*1e-5;
        return findScalarRootBrent(x_[lookUpIndex-1]-eps, x_[lookUpIndex]+eps, localPolynomial);
    }

    auto lookUpIndex_(const std::vector<Scalar>& vec, const Scalar v, bool increasing) const
    {
        return increasing ? lookUpIndexIncreasing_(vec, v) : lookUpIndexDecreasing_(vec, v);
    }

    auto lookUpIndexIncreasing_(const std::vector<Scalar>& vec, const Scalar v) const
    {
        const auto lookUpIndex = std::distance(vec.begin(), std::lower_bound(vec.begin(), vec.end(), v));
        assert(lookUpIndex != 0 && lookUpIndex < vec.size());
        return lookUpIndex;
    }

    auto lookUpIndexDecreasing_(const std::vector<Scalar>& vec, const Scalar v) const
    {
        const auto lookUpIndex = vec.size() - std::distance(vec.rbegin(), std::upper_bound(vec.rbegin(), vec.rend(), v));
        assert(lookUpIndex != 0 && lookUpIndex < vec.size());
        return lookUpIndex;
    }

    std::vector<Scalar> x_; //!< the x-coordinates
    std::vector<Scalar> y_; //!< the y-coordinates
    std::vector<Scalar> m_; //!< the slope for each control point
    std::size_t numPoints_; //!< the number of control points
    bool increasingX_; //!< if we are increasing monotone or not
    bool increasingY_; //!< if we are increasing monotone or not
};

} // end namespace Dumux

#endif
