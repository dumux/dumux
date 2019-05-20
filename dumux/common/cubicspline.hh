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
 * \brief A simple implementation of a cubic spline
 */
#ifndef DUMUX_COMMON_CUBIC_SPLINE_HH
#define DUMUX_COMMON_CUBIC_SPLINE_HH

#include <iterator>
#include <vector>
#include <algorithm>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/istl/btdmatrix.hh>
#include <dune/istl/bvector.hh>

namespace Dumux {

/*!
 * \ingroup Common
 * \brief A simple implementation of a natural cubic spline
 * \note We follow the notation at http://mathworld.wolfram.com/CubicSpline.html
 */
template<class Scalar = double>
class CubicSpline
{
public:
    /*!
     * \brief Default constructor
     */
    CubicSpline() = default;

    /*!
     * \brief Contruct a natural cubic spline from the control points (x[i], y[i])
     * \param x a vector of x-coordinates
     * \param y a vector of y-coordinates
     */
    CubicSpline(const std::vector<Scalar>& x, const std::vector<Scalar> y)
    {
        // check some requirements
        assert (x.size() == y.size());
        assert (x.size() >=2);
        assert (std::is_sorted(x.begin(), x.end()));

        updatePoints(x, y);
    }

    /*!
     * \brief Create a natural cubic spline from the control points (x[i], y[i])
     * \note we enforce continuous second derivatives in the inside and zero second derivatives at the boundary
     * \param x a vector of x-coordinates
     * \param y a vector of y-coordinates
     */
    void updatePoints(const std::vector<Scalar>& x, const std::vector<Scalar>& y)
    {
        // save a copy of the control points
        x_ = x;

        // the number of control points
        numPoints_ = x.size();

        const auto numSegments = numPoints_-1;
        coeff_.resize(numSegments*4+2); // 4 coefficients for each segment + last y-value and derivative

        // construct a block-tridiagonal matrix system solving for the derivatives
        Dune::BTDMatrix<Dune::FieldMatrix<double, 1, 1>> matrix(numPoints_);
        Dune::BlockVector<Dune::FieldVector<double, 1>> rhs(numPoints_);
        Dune::BlockVector<Dune::FieldVector<double, 1>> d(numPoints_);

        // assemble matrix and rhs row-wise
        matrix[0][0] = 2.0;
        matrix[0][1] = 1.0;
        rhs[0] = 3.0*(y[1] - y[0]);
        for (int i = 1; i < numPoints_-1; ++i)
        {
            matrix[i][i-1] = 1.0;
            matrix[i][i] = 4.0;
            matrix[i][i+1] = 1.0;
            rhs[i] = 3.0*(y[i+1] - y[i-1]);
        }
        matrix[numPoints_-1][numPoints_-1] = 2.0;
        matrix[numPoints_-1][numPoints_-2] = 1.0;
        rhs[numPoints_-1] = 3.0*(y[numPoints_-1] - y[numPoints_-2]);

        // solve for derivatives
        matrix.solve(d, rhs);

        // compute coefficients
        std::size_t offset = 0;
        for (int i = 0; i < numSegments; ++i, offset += 4)
        {
            coeff_[offset+0] = y[i];
            coeff_[offset+1] = d[i];
            coeff_[offset+2] = 3.0*(y[i+1]-y[i]) - 2.0*d[i] - d[i+1];
            coeff_[offset+3] = 2.0*(y[i]-y[i+1]) + d[i] + d[i+1];
        }
        coeff_[offset+0] = y[numPoints_-1];
        coeff_[offset+1] = d[numPoints_-1];
    }

    /*!
     * \brief Evaluate the y value at a given x value
     * \param x the x-coordinate
     * \note We extrapolate linearly if out of bounds
     */
    Scalar eval(const Scalar x) const
    {
        if (x <= x_[0])
            return coeff_[0] + evalDerivative(x_[0])*(x - x_[0]);
        else if (x > x_[numPoints_-1])
            return coeff_[(numPoints_-1)*4] + evalDerivative(x_[numPoints_-1])*(x - x_[numPoints_-1]);
        else
        {
            const auto lookUpIndex = std::distance(x_.begin(), std::lower_bound(x_.begin(), x_.end(), x));
            assert(lookUpIndex != 0);

            // get coefficients
            const auto* coeff = coeff_.data() + (lookUpIndex-1)*4;

            // interpolate parametrization parameter t in [0,1]
            const auto t = (x - x_[lookUpIndex-1])/(x_[lookUpIndex] - x_[lookUpIndex-1]);
            return coeff[0] + t*(coeff[1] + t*coeff[2] + t*t*coeff[3]);
        }
    }

    /*!
     * \brief Evaluate the first derivative dy/dx at a given x value
     * \param x the x-coordinate
     * \note We extrapolate linearly if out of bounds
     */
    Scalar evalDerivative(const Scalar x) const
    {
        if (x <= x_[0])
            return coeff_[1]/(x_[1] - x_[0]);
        else if (x > x_[numPoints_-1])
            return coeff_[(numPoints_-1)*4 + 1]/(x_[numPoints_-1] - x_[numPoints_-2]);
        else
        {
            const auto lookUpIndex = std::distance(x_.begin(), std::lower_bound(x_.begin(), x_.end(), x));
            assert(lookUpIndex != 0);

            // get coefficients
            const auto* coeff = coeff_.data() + (lookUpIndex-1)*4;

            // interpolate parametrization parameter t in [0,1]
            const auto t = (x - x_[lookUpIndex-1])/(x_[lookUpIndex] - x_[lookUpIndex-1]);
            const auto dtdx = 1.0/(x_[lookUpIndex] - x_[lookUpIndex-1]);
            return dtdx*(coeff[1] + t*(2.0*coeff[2] + t*3.0*coeff[3]));
        }
    }

private:
    std::vector<Scalar> x_; //!< the x-coordinates
    std::size_t numPoints_; //!< the number of control points
    std::vector<Scalar> coeff_; //!< the spline coefficients
};

} // end namespace Dumux

#endif
