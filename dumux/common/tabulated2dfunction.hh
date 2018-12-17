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
 * \brief Implements tabulation for a two-dimensional function.
 */
#ifndef DUMUX_TABULATED_2D_FUNCTION_HH
#define DUMUX_TABULATED_2D_FUNCTION_HH

#include <vector>
#include <cassert>

namespace Dumux {

/*!
 * \ingroup Common
 * \brief Implements tabulation for a two-dimensional function.
 *
 * This class can be used to tabulate a two dimensional function
 * \f$f(x, y)\f$ over the range \f$[x_{min}, x_{max}] \times [y_{min},
 * y_{max}]\f$. For this, the ranges of the \f$x\f$ and \f$y\f$ axes are
 * divided into \f$m\f$ and \f$n\f$ sub-intervals and the values of
 * \f$f(x_i, y_j)\f$ need to be provided. Here, \f$x_i\f$ and
 * \f$y_j\f$ are the largest positions of the \f$i\f$-th and
 * \f$j\f$-th intervall. Between these sampling points this tabulation
 * class uses linear interpolation.
 */
template<class Scalar>
class Tabulated2DFunction
{
public:
    /*!
     * \brief Default constructor.
     */
    Tabulated2DFunction()
    { };

    /*!
     * \brief Constructor where the tabulation parameters are already
     *        provided.
     */
    Tabulated2DFunction(Scalar xMin, Scalar xMax, int m,
                        Scalar yMin, Scalar yMax, int n)
    {
        resize(xMin, xMax, m, yMin, yMax, n);
    };

    /*!
     * \brief Resize the tabulation to a new range.
     *
     * \note _All_ previously specified sampling points become invalid
     *       after calling this method. You have to supply new ones.
     */
    void resize(Scalar xMin, Scalar xMax, int m,
                Scalar yMin, Scalar yMax, int n)
    {
        samples_.resize(m*n);

        m_ = m;
        n_ = n;

        xMin_ = xMin;
        xMax_ = xMax;

        yMin_ = yMin;
        yMax_ = yMax;
    }

    /*!
     * \brief Return the position on the x-axis of the i-th interval.
     */
    Scalar iToX(int i) const
    {
        assert(0 <= i && i < m_);

        return xMin_ + i*(xMax_ - xMin_)/(m_ - 1);
    }

    /*!
     * \brief Return the position on the y-axis of the j-th interval.
     */
    Scalar jToY(int j) const
    {
        assert(0 <= j && j < n_);

        return yMin_ + j*(yMax_ - yMin_)/(n_ - 1);
    }

    /*!
     * \brief Return the interval index of a given position on the x-axis.
     *
     * This method returns a *floating point* number. The integer part
     * should be interpreted as intervall, the decimal places are the
     * position of the x value between the i-th and the (i+1)-th
     * sample point.
     */
    Scalar xToI(Scalar x) const
    {
        return (x - xMin_)/(xMax_ - xMin_)*m_;
    }


    /*!
     * \brief Return the interval index of a given position on the y-axis.
     *
     * This method returns a *floating point* number. The integer part
     * should be interpreted as intervall, the decimal places are the
     * position of the y value between the j-th and the (j+1)-th
     * sample point.
     */
    Scalar yToJ(Scalar y) const
    {
        return (y - yMin_)/(yMax_ - yMin_)*n_;
    }


    /*!
     * \brief Get the value of the sample point which is at the
     *         intersection of the \f$i\f$-th interval of the x-Axis
     *         and the \f$j\f$-th of the y-Axis.
     */
    Scalar getSamplePoint(int i, int j) const
    {
        assert(0 <= i && i < m_);
        assert(0 <= j && j < n_);

        return samples_[j*m_ + i];
    }

    /*!
     * \brief Set the value of the sample point which is at the
     *        intersection of the \f$i\f$-th interval of the x-Axis
     *        and the \f$j\f$-th of the y-Axis.
     */
    void setSamplePoint(int i, int j, Scalar value)
    {
        assert(0 <= i && i < m_);
        assert(0 <= j && j < n_);

        samples_[j*m_ + i] = value;
    }

    /*!
     * \brief Return an interpolated value.
     */
    Scalar get(Scalar x, Scalar y) const
    {
        Scalar alpha = xToI(x);
        Scalar beta = yToJ(y);

        using std::max;
        using std::min;
        int i = max(0, min(m_, static_cast<int>(alpha)));
        int j = max(0, min(n_, static_cast<int>(beta)));

        alpha -= i;
        beta -= j;

        // bi-linear interpolation
        Scalar s1 = getSamplePoint(i, j)*(1.0 - alpha) + getSamplePoint(i + 1, j)*alpha;
        Scalar s2 = getSamplePoint(i, j + 1)*(1.0 - alpha) + getSamplePoint(i + 1, j + 1)*alpha;
        return s1*(1.0 - beta) + s2*beta;
    }

    /*!
     * \brief The () operator
     *
     * This is just a convenience alias for get(x,y);
     */
    Scalar operator()(Scalar x, Scalar y) const
    { return get(x, y); }

private:
    // the vector which contains the values of the sample points
    // f(x_i, y_j). don't use this directly, use getSamplePoint(i,j)
    // instead!
    std::vector<Scalar> samples_;

    // the number of sample points in x direction
    int m_;

    // the number of sample points in y direction
    int n_;

    // the range of the tabulation on the x axis
    Scalar xMin_;
    Scalar xMax_;

    // the range of the tabulation on the y axis
    Scalar yMin_;
    Scalar yMax_;

};

} // end namespace Dumux

#endif
