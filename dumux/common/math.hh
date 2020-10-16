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
 * \brief Define some often used mathematical functions
 */
#ifndef DUMUX_MATH_HH
#define DUMUX_MATH_HH

#include <algorithm>
#include <cmath>
#include <utility>

#include <dune/common/typetraits.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/common/float_cmp.hh>

namespace Dumux {

/*!
 * \ingroup Common
 * \brief Calculate the (weighted) arithmetic mean of two scalar values.
 *
 * \param x The first input value
 * \param y The second input value
 * \param wx The first weight
 * \param wy The second weight
 */
template <class Scalar>
constexpr Scalar arithmeticMean(Scalar x, Scalar y, Scalar wx = 1.0, Scalar wy = 1.0) noexcept
{
    static_assert(Dune::IsNumber<Scalar>::value, "The arguments x, y, wx, and wy have to be numbers!");

    if (x*y <= 0)
        return 0;
    return (x * wx + y * wy)/(wx + wy);
}

/*!
 * \ingroup Common
 * \brief Calculate the (weighted) harmonic mean of two scalar values.
 *
 * \param x The first input value
 * \param y The second input value
 * \param wx The first weight
 * \param wy The second weight
 */
template <class Scalar>
constexpr Scalar harmonicMean(Scalar x, Scalar y, Scalar wx = 1.0, Scalar wy = 1.0) noexcept
{
    static_assert(Dune::IsNumber<Scalar>::value, "The arguments x, y, wx, and wy have to be numbers!");

    if (x*y <= 0)
        return 0;
    return (wx + wy) * x * y / (wy * x + wx * y);
}

/*!
 * \ingroup Common
 * \brief Calculate the geometric mean of two scalar values.
 *
 * \param x The first input value
 * \param y The second input value
 * \note as std::sqrt is not constexpr this function is not constexpr
 */
template <class Scalar>
Scalar geometricMean(Scalar x, Scalar y) noexcept
{
    static_assert(Dune::IsNumber<Scalar>::value, "The arguments x and y have to be numbers!");

    if (x*y <= 0)
        return 0;
    using std::sqrt;
    return sqrt(x*y)*sign(x);
}

/*!
 * \ingroup Common
 * \brief Calculate the harmonic mean of a fixed-size matrix.
 *
 * This is done by calculating the harmonic mean for each entry
 * individually. The result is stored the first argument.
 *
 * \param K The matrix for storing the result
 * \param Ki The first input matrix
 * \param Kj The second input matrix
 */
template <class Scalar, int m, int n>
void harmonicMeanMatrix(Dune::FieldMatrix<Scalar, m, n> &K,
                        const Dune::FieldMatrix<Scalar, m, n> &Ki,
                        const Dune::FieldMatrix<Scalar, m, n> &Kj)
{
    for (int rowIdx=0; rowIdx < m; rowIdx++){
        for (int colIdx=0; colIdx< n; colIdx++){
            if (Dune::FloatCmp::ne<Scalar>(Ki[rowIdx][colIdx], Kj[rowIdx][colIdx])) {
                K[rowIdx][colIdx] =
                    harmonicMean(Ki[rowIdx][colIdx],
                                 Kj[rowIdx][colIdx]);
            }
            else
                K[rowIdx][colIdx] = Ki[rowIdx][colIdx];
        }
    }
}

/*!
 * \ingroup Common
 * \brief Invert a linear polynomial analytically
 *
 * The polynomial is defined as
 * \f[ p(x) = a\; x + b \f]
 *
 * This method Returns the number of solutions which are in the real
 * numbers, i.e. 1 except if the slope of the line is 0.
 *
 * \param sol Container into which the solutions are written
 * \param a The coefficiont for the linear term
 * \param b The coefficiont for the constant term
 */
template <class Scalar, class SolContainer>
int invertLinearPolynomial(SolContainer &sol,
                           Scalar a,
                           Scalar b)
{
    if (a == 0.0)
        return 0;

    sol[0] = -b/a;
    return 1;
}

/*!
 * \ingroup Common
 * \brief Invert a quadratic polynomial analytically
 *
 * The polynomial is defined as
 * \f[ p(x) = a\; x^2 + + b\;x + c \f]
 *
 * This method returns the number of solutions which are in the real
 * numbers. The "sol" argument contains the real roots of the parabola
 * in order with the smallest root first.
 *
 * \param sol Container into which the solutions are written
 * \param a The coefficiont for the quadratic term
 * \param b The coefficiont for the linear term
 * \param c The coefficiont for the constant term
 */
template <class Scalar, class SolContainer>
int invertQuadraticPolynomial(SolContainer &sol,
                              Scalar a,
                              Scalar b,
                              Scalar c)
{
    // check for a line
    if (a == 0.0)
        return invertLinearPolynomial(sol, b, c);

    // discriminant
    Scalar Delta = b*b - 4*a*c;
    if (Delta < 0)
        return 0; // no real roots

    using std::sqrt;
    Delta = sqrt(Delta);
    sol[0] = (- b + Delta)/(2*a);
    sol[1] = (- b - Delta)/(2*a);

    // sort the result
    if (sol[0] > sol[1])
    {
        using std::swap;
        swap(sol[0], sol[1]);
    }
    return 2; // two real roots
}

//! \cond false
template <class Scalar, class SolContainer>
void invertCubicPolynomialPostProcess_(SolContainer &sol,
                                       int numSol,
                                       Scalar a,
                                       Scalar b,
                                       Scalar c,
                                       Scalar d)
{
    // do one Newton iteration on the analytic solution if the
    // precision is increased
    for (int i = 0; i < numSol; ++i) {
        Scalar x = sol[i];
        Scalar fOld = d + x*(c + x*(b + x*a));

        Scalar fPrime = c + x*(2*b + x*3*a);
        if (fPrime == 0.0)
            continue;
        x -= fOld/fPrime;

        Scalar fNew = d + x*(c + x*(b + x*a));
        using std::abs;
        if (abs(fNew) < abs(fOld))
            sol[i] = x;
    }
}
//! \endcond

/*!
 * \ingroup Common
 * \brief Invert a cubic polynomial analytically
 *
 * The polynomial is defined as
 * \f[ p(x) = a\; x^3 + + b\;x^3 + c\;x + d \f]
 *
 * This method teturns the number of solutions which are in the real
 * numbers. The "sol" argument contains the real roots of the cubic
 * polynomial in order with the smallest root first.
 *
 * \param sol Container into which the solutions are written
 * \param a The coefficiont for the cubic term
 * \param b The coefficiont for the quadratic term
 * \param c The coefficiont for the linear term
 * \param d The coefficiont for the constant term
 */
template <class Scalar, class SolContainer>
int invertCubicPolynomial(SolContainer *sol,
                          Scalar a,
                          Scalar b,
                          Scalar c,
                          Scalar d)
{
    // reduces to a quadratic polynomial
    if (a == 0)
        return invertQuadraticPolynomial(sol, b, c, d);

    // normalize the polynomial
    b /= a;
    c /= a;
    d /= a;
    a = 1;

    // get rid of the quadratic term by subsituting x = t - b/3
    Scalar p = c - b*b/3;
    Scalar q = d + (2*b*b*b - 9*b*c)/27;

    // now we are at the form t^3 + p*t + q = 0. First we handle some
    // special cases to avoid divisions by zero later...
    if (p == 0.0 && q == 0.0) {
        // t^3 = 0, i.e. triple root at t = 0
        sol[0] = sol[1] = sol[2] = 0.0 - b/3;
        return 3;
    }
    else if (p == 0.0 && q != 0.0) {
        // t^3 + q = 0,
        //
        // i. e. single real root at t=curt(q)
        using std::cbrt;
        Scalar t = cbrt(q);
        sol[0] = t - b/3;

        return 1;
    }
    else if (p != 0.0 && q == 0.0) {
        // t^3 + p*t = 0 = t*(t^2 + p),
        //
        // i. e. roots at t = 0, t^2 + p = 0
        if (p > 0) {
            sol[0] = 0.0 - b/3;
            return 1; // only a single real root at t=0
        }

        // two additional real roots at t = sqrt(-p) and t = -sqrt(-p)
        using std::sqrt;
        sol[0] = -sqrt(-p) - b/3;
        sol[1] = 0.0 - b/3;
        sol[2] = sqrt(-p) - b/3;

        return 3;
    }

    // At this point
    //
    // t^3 + p*t + q = 0
    //
    // with p != 0 and q != 0 holds. Introducing the variables u and v
    // with the properties
    //
    //   u + v = t       and       3*u*v + p = 0
    //
    // leads to
    //
    // u^3 + v^3 + q = 0 .
    //
    // multiplying both sides with u^3 and taking advantage of the
    // fact that u*v = -p/3 leads to
    //
    // u^6 + q*u^3 - p^3/27 = 0
    //
    // Now, substituting u^3 = w yields
    //
    // w^2 + q*w - p^3/27 = 0
    //
    // This is a quadratic equation with the solutions
    //
    // w = -q/2 + sqrt(q^2/4 + p^3/27) and
    // w = -q/2 - sqrt(q^2/4 + p^3/27)
    //
    // Since w is equivalent to u^3 it is sufficient to only look at
    // one of the two cases. Then, there are still 2 cases: positive
    // and negative discriminant.
    Scalar wDisc = q*q/4 + p*p*p/27;
    if (wDisc >= 0) { // the positive discriminant case:
        // calculate the cube root of - q/2 + sqrt(q^2/4 + p^3/27)
        using std::cbrt;
        using std::sqrt;
        Scalar u = cbrt(-q/2 + sqrt(wDisc));

        // at this point, u != 0 since p^3 = 0 is necessary in order
        // for u = 0 to hold, so
        sol[0] = u - p/(3*u) - b/3;
        // does not produce a division by zero. the remaining two
        // roots of u are rotated by +- 2/3*pi in the complex plane
        // and thus not considered here
        invertCubicPolynomialPostProcess_(sol, 1, a, b, c, d);
        return 1;
    }
    else { // the negative discriminant case:
        Scalar uCubedRe = - q/2;
        using std::sqrt;
        Scalar uCubedIm = sqrt(-wDisc);
        // calculate the cube root of - q/2 + sqrt(q^2/4 + p^3/27)
        using std::cbrt;
        Scalar uAbs  = cbrt(sqrt(uCubedRe*uCubedRe + uCubedIm*uCubedIm));
        using std::atan2;
        Scalar phi = atan2(uCubedIm, uCubedRe)/3;

        // calculate the length and the angle of the primitive root

        // with the definitions from above it follows that
        //
        // x = u - p/(3*u) - b/3
        //
        // where x and u are complex numbers. Rewritten in polar form
        // this is equivalent to
        //
        // x = |u|*e^(i*phi) - p*e^(-i*phi)/(3*|u|) - b/3 .
        //
        // Factoring out the e^ terms and subtracting the additional
        // terms, yields
        //
        // x = (e^(i*phi) + e^(-i*phi))*(|u| - p/(3*|u|)) - y - b/3
        //
        // with
        //
        // y = - |u|*e^(-i*phi) + p*e^(i*phi)/(3*|u|) .
        //
        // The crucial observation is the fact that y is the conjugate
        // of - x + b/3. This means that after taking advantage of the
        // relation
        //
        // e^(i*phi) + e^(-i*phi) = 2*cos(phi)
        //
        // the equation
        //
        // x = 2*cos(phi)*(|u| - p / (3*|u|)) - conj(x) - 2*b/3
        //
        // holds. Since |u|, p, b and cos(phi) are real numbers, it
        // follows that Im(x) = - Im(x) and thus Im(x) = 0. This
        // implies
        //
        // Re(x) = x = cos(phi)*(|u| - p / (3*|u|)) - b/3 .
        //
        // Considering the fact that u is a cubic root, we have three
        // values for phi which differ by 2/3*pi. This allows to
        // calculate the three real roots of the polynomial:
        for (int i = 0; i < 3; ++i) {
            using std::cos;
            sol[i] = cos(phi)*(uAbs - p/(3*uAbs)) - b/3;
            phi += 2*M_PI/3;
        }

        // post process the obtained solution to increase numerical
        // precision
        invertCubicPolynomialPostProcess_(sol, 3, a, b, c, d);

        // sort the result
        using std::sort;
        sort(sol, sol + 3);

        return 3;
    }

    // NOT REACHABLE!
    return 0;
}

/*!
 * \ingroup Common
 * \brief Comparison of two position vectors
 *
 * Compares an current position vector with a reference vector, and returns true
 * if the position vector is larger.
 * "Larger" in this case means that all the entries of each spacial dimension are
 * larger compared to the reference vector.
 *
 * \param pos Vektor holding the current Position that is to be checked
 * \param smallerVec Reference vector, holding the minimum values for comparison.
 */
template <class Scalar, int dim>
bool isLarger(const Dune::FieldVector<Scalar, dim> &pos,
              const Dune::FieldVector<Scalar, dim> &smallerVec)
{
    for (int i=0; i < dim; i++)
    {
        if (pos[i]<= smallerVec[i])
        {
            return false;
        }
    }
    return true;
}

/*!
 * \ingroup Common
 * \brief Comparison of two position vectors
 *
 * Compares an current position vector with a reference vector, and returns true
 * if the position vector is smaller.
 * "Smaller" in this case means that all the entries of each spacial dimension are
 * smaller in comparison with the reference vector.
 *
 * \param pos Vektor holding the current Position that is to be checked
 * \param largerVec Reference vector, holding the maximum values for comparison.
 */
template <class Scalar, int dim>
bool isSmaller(const Dune::FieldVector<Scalar, dim> &pos,
               const Dune::FieldVector<Scalar, dim> &largerVec)
{
    for (int i=0; i < dim; i++)
    {
        if (pos[i]>= largerVec[i])
        {
            return false;
        }
    }
    return true;
}

/*!
 * \ingroup Common
 * \brief Comparison of three position vectors
 *
 * Compares an current position vector with two reference vector, and returns true
 * if the position vector lies in between them.
 * "Between" in this case means that all the entries of each spacial dimension are
 * smaller in comparison with the larger reference vector as well as larger campared
 * to the smaller reference.
 * This is comfortable to cheack weather the current position is located inside or
 * outside of a lense with different properties.
 *
 * \param pos Vektor holding the current Position that is to be checked
 * \param smallerVec Reference vector, holding the minimum values for comparison.
 * \param largerVec Reference vector, holding the maximum values for comparison.
 */
template <class Scalar, int dim>
bool isBetween(const Dune::FieldVector<Scalar, dim> &pos,
              const Dune::FieldVector<Scalar, dim> &smallerVec,
              const Dune::FieldVector<Scalar, dim> &largerVec)
{
   if (isLarger(pos, smallerVec) && isSmaller(pos, largerVec))
       {
           return true;
       }
   else
       return false;
}

//! forward declaration of the linear interpolation policy (default)
namespace InterpolationPolicy { struct Linear; }

/*!
 * \ingroup Common
 * \brief a generic function to interpolate given a set of parameters and an interpolation point
 * \param params the parameters used for interpolation (depends on the policy used)
 * \param ip the interpolation point
 */
template <class Policy = InterpolationPolicy::Linear, class Scalar, class ... Parameter>
Scalar interpolate(Scalar ip, Parameter&& ... params)
{ return Policy::interpolate(ip, std::forward<Parameter>(params) ...); }

/*!
 * \ingroup Common
 * \brief Interpolation policies
 */
namespace InterpolationPolicy  {

/*!
 * \ingroup Common
 * \brief interpolate linearly between two given values
 */
struct Linear
{
    /*!
     * \brief interpolate linearly between two given values
     * \param ip the interpolation point in [0,1]
     * \param params array with the lower and upper bound
     */
    template<class Scalar>
    static constexpr Scalar interpolate(Scalar ip, const std::array<Scalar, 2>& params)
    {
        return params[0]*(1.0 - ip) + params[1]*ip;
    }
};

/*!
 * \ingroup Common
 * \brief interpolate linearly in a piecewise linear function (tabularized function)
 */
struct LinearTable
{
    /*!
     * \brief interpolate linearly in a piecewise linear function (tabularized function)
     * \param ip the interpolation point
     * \param range positions of values
     * \param values values to interpolate from
     * \note if the interpolation point is out of bounds this will return the bounds
     */
    template<class Scalar, class RandomAccessContainer0, class RandomAccessContainer1>
    static constexpr Scalar interpolate(Scalar ip, const RandomAccessContainer0& range, const RandomAccessContainer1& values)
    {
        // check bounds
        if (ip > range.back()) return values.back();
        if (ip < range[0]) return values[0];

        // if we are within bounds find the index of the lower bound
        const auto lookUpIndex = std::distance(range.begin(), std::lower_bound(range.begin(), range.end(), ip));
        if (lookUpIndex == 0)
            return values[0];

        const auto ipLinear = (ip - range[lookUpIndex-1])/(range[lookUpIndex] - range[lookUpIndex-1]);
        return Dumux::interpolate<Linear>(ipLinear, std::array<Scalar, 2>{{values[lookUpIndex-1], values[lookUpIndex]}});
    }

    template<class Scalar, class RandomAccessContainer>
    static constexpr Scalar interpolate(Scalar ip, const std::pair<RandomAccessContainer, RandomAccessContainer>& table)
    {
        const auto& [range, values] = table;
        return interpolate(ip, range, values);
    }
};

} // end namespace InterpolationPolicy

/*!
 * \ingroup Common
 * \brief Generates linearly spaced vectors
 *
 * \param begin The first value in the vector
 * \param end The last value in the vector
 * \param samples The size of the vector
 * \param endPoint if the range is including the interval's end point or not
 */
template <class Scalar>
std::vector<Scalar> linspace(const Scalar begin, const Scalar end,
                             std::size_t samples,
                             bool endPoint = true)
{
    using std::max;
    samples = max(std::size_t{2}, samples); // only makes sense for 2 or more samples
    const Scalar divisor = endPoint ? samples-1 : samples;
    const Scalar delta = (end-begin)/divisor;
    std::vector<Scalar> vec(samples);
    for (std::size_t i = 0; i < samples; ++i)
        vec[i] = begin + i*delta;
    return vec;
}


/*!
 * \ingroup Common
 * \brief Evaluates the Antoine equation used to calculate the vapour
 *        pressure of various liquids.
 *
 * See http://en.wikipedia.org/wiki/Antoine_equation
 *
 * \param temperature The temperature [K] of the fluid
 * \param A The first coefficient for the Antoine equation
 * \param B The first coefficient for the Antoine equation
 * \param C The first coefficient for the Antoine equation
 */
template <class Scalar>
Scalar antoine(Scalar temperature,
               Scalar A,
               Scalar B,
               Scalar C)
{
    const Scalar ln10 = 2.3025850929940459;
    using std::exp;
    return exp(ln10*(A - B/(C + temperature)));
}

/*!
 * \ingroup Common
 * \brief Sign or signum function.
 *
 * Returns 1 for a positive argument.
 * Returns -1 for a negative argument.
 * Returns 0 if the argument is zero.
 */
template<class ValueType>
constexpr int sign(const ValueType& value) noexcept
{
    return (ValueType(0) < value) - (value < ValueType(0));
}

/*!
 * \ingroup Common
 * \brief Cross product of two vectors in three-dimensional Euclidean space
 *
 * \param vec1 The first vector
 * \param vec2 The second vector
 */
template <class Scalar>
Dune::FieldVector<Scalar, 3> crossProduct(const Dune::FieldVector<Scalar, 3> &vec1,
                                          const Dune::FieldVector<Scalar, 3> &vec2)
{
    return {vec1[1]*vec2[2]-vec1[2]*vec2[1],
            vec1[2]*vec2[0]-vec1[0]*vec2[2],
            vec1[0]*vec2[1]-vec1[1]*vec2[0]};
}

/*!
 * \ingroup Common
 * \brief Cross product of two vectors in two-dimensional Euclidean space retuning scalar
 *
 * \param vec1 The first vector
 * \param vec2 The second vector
 */
template <class Scalar>
Scalar crossProduct(const Dune::FieldVector<Scalar, 2> &vec1,
                    const Dune::FieldVector<Scalar, 2> &vec2)
{   return vec1[0]*vec2[1]-vec1[1]*vec2[0]; }

/*!
 * \ingroup Common
 * \brief Triple product of three vectors in three-dimensional Euclidean space retuning scalar
 *
 * \param vec1 The first vector
 * \param vec2 The second vector
 * \param vec3 The third vector
 */
template <class Scalar>
Scalar tripleProduct(const Dune::FieldVector<Scalar, 3> &vec1,
                     const Dune::FieldVector<Scalar, 3> &vec2,
                     const Dune::FieldVector<Scalar, 3> &vec3)
{   return crossProduct<Scalar>(vec1, vec2)*vec3; }

/*!
 * \ingroup Common
 * \brief Transpose a FieldMatrix
 *
 * \param M The matrix to be transposed
 */
template <class Scalar, int m, int n>
Dune::FieldMatrix<Scalar, n, m> getTransposed(const Dune::FieldMatrix<Scalar, m, n>& M)
{
    Dune::FieldMatrix<Scalar, n, m> T;
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < n; ++j)
            T[j][i] = M[i][j];

    return T;
}

/*!
 * \ingroup Common
 * \brief Transpose a DynamicMatrix
 *
 * \param M The matrix to be transposed
 */
template <class Scalar>
Dune::DynamicMatrix<Scalar> getTransposed(const Dune::DynamicMatrix<Scalar>& M)
{
    std::size_t rows_T = M.M();
    std::size_t cols_T = M.N();

    Dune::DynamicMatrix<Scalar> M_T(rows_T, cols_T, 0.0);

    for (std::size_t i = 0; i < rows_T; ++i)
        for (std::size_t j = 0; j < cols_T; ++j)
            M_T[i][j] = M[j][i];

    return M_T;
}

/*!
 * \ingroup Common
 * \brief Multiply two dynamic matrices
 *
 * \param M1 The first dynamic matrix
 * \param M2 The second dynamic matrix (to be multiplied to M1 from the right side)
 */
template <class Scalar>
Dune::DynamicMatrix<Scalar> multiplyMatrices(const Dune::DynamicMatrix<Scalar> &M1,
                                             const Dune::DynamicMatrix<Scalar> &M2)
{
    using size_type = typename Dune::DynamicMatrix<Scalar>::size_type;
    const size_type rows = M1.N();
    const size_type cols = M2.M();

    DUNE_ASSERT_BOUNDS(M1.M() == M2.N());

    Dune::DynamicMatrix<Scalar> result(rows, cols, 0.0);
    for (size_type i = 0; i < rows; i++)
        for (size_type j = 0; j < cols; j++)
            for (size_type k = 0; k < M1.M(); k++)
                result[i][j] += M1[i][k]*M2[k][j];

    return result;
}

/*!
 * \ingroup Common
 * \brief Multiply two field matrices
 *
 * \param M1 The first field matrix
 * \param M2 The second field matrix (to be multiplied to M1 from the right side)
 */
template <class Scalar, int rows1, int cols1, int cols2>
Dune::FieldMatrix<Scalar, rows1, cols2> multiplyMatrices(const Dune::FieldMatrix<Scalar, rows1, cols1> &M1,
                                                         const Dune::FieldMatrix<Scalar, cols1, cols2> &M2)
{
    using size_type = typename Dune::FieldMatrix<Scalar, rows1, cols2>::size_type;

    Dune::FieldMatrix<Scalar, rows1, cols2> result(0.0);
    for (size_type i = 0; i < rows1; i++)
        for (size_type j = 0; j < cols2; j++)
            for (size_type k = 0; k < cols1; k++)
                result[i][j] += M1[i][k]*M2[k][j];

    return result;
}


/*!
 * \ingroup Common
 * \brief Trace of a dense matrix
 *
 * \param M The dense matrix
 */
template <class MatrixType>
typename Dune::DenseMatrix<MatrixType>::field_type
trace(const Dune::DenseMatrix<MatrixType>& M)
{
    const auto rows = M.N();
    DUNE_ASSERT_BOUNDS(rows == M.M()); // rows == cols

    using MatType = Dune::DenseMatrix<MatrixType>;
    typename MatType::field_type trace = 0.0;

    for (typename MatType::size_type i = 0; i < rows; ++i)
        trace += M[i][i];

    return trace;
}

/*!
 * \ingroup Common
 * \brief Returns the result of the projection of
 *        a vector v with a Matrix M.
 *
 *        Note: We use DenseVector and DenseMatrix here so that
 *        it can be used with the statically and dynamically
 *        allocated Dune Vectors/Matrices. Size mismatch
 *        assertions are done in the respective Dune classes.
 *
 * \param M The matrix
 * \param v The vector
 */
template <class MAT, class V>
typename Dune::DenseVector<V>::derived_type
mv(const Dune::DenseMatrix<MAT>& M,
   const Dune::DenseVector<V>& v)
{
    typename Dune::DenseVector<V>::derived_type res(v);
    M.mv(v, res);
    return res;
}

/*!
 * \ingroup Common
 * \brief Returns the result of a vector v multiplied by a scalar m.
 *
 *        Note: We use DenseVector and DenseMatrix here so that
 *        it can be used with the statically and dynamically
 *        allocated Dune Vectors/Matrices. Size mismatch
 *        assertions are done in the respective Dune classes.
 *
 * \note We need the enable_if to make sure that only Scalars
 *       fit here. Matrix types are forced to use the above
 *       mv(DenseMatrix, DenseVector) instead.
 *
 * \param m The scale factor
 * \param v The vector
 */

template <class FieldScalar, class V>
typename std::enable_if_t<Dune::IsNumber<FieldScalar>::value,
                          typename Dune::DenseVector<V>::derived_type>
mv(const FieldScalar m, const Dune::DenseVector<V>& v)
{
    typename Dune::DenseVector<V>::derived_type res(v);
    res *= m;
    return res;
}

/*!
 * \ingroup Common
 * \brief Evaluates the scalar product of a vector v2, projected by
 *        a matrix M, with a vector v1.
 *
 *        Note: We use DenseVector and DenseMatrix here so that
 *        it can be used with the statically and dynamically
 *        allocated Dune Vectors/Matrices. Size mismatch
 *        assertions are done in the respective Dune classes.
 *
 * \param v1 The first vector
 * \param M The matrix
 * \param v2 The second vector
 */
template <class V1, class MAT, class V2>
typename Dune::DenseMatrix<MAT>::value_type
vtmv(const Dune::DenseVector<V1>& v1,
     const Dune::DenseMatrix<MAT>& M,
     const Dune::DenseVector<V2>& v2)
{
    return v1*mv(M, v2);
}

/*!
 * \ingroup Common
 * \brief Evaluates the scalar product of a vector v2, scaled by
 *        a scalar m, with a vector v1.
 *
 *        Note: We use DenseVector and DenseMatrix here so that
 *        it can be used with the statically and dynamically
 *        allocated Dune Vectors/Matrices. Size mismatch
 *        assertions are done in the respective Dune classes.
 *
 * \note We need the enable_if to make sure that only Scalars
 *       fit here. Matrix types are forced to use the above
 *       vtmv(DenseVector, DenseMatrix, DenseVector) instead.
 *
 * \param v1 The first vector
 * \param m The scale factor
 * \param v2 The second vector
 */

template <class V1, class FieldScalar, class V2>
typename std::enable_if_t<Dune::IsNumber<FieldScalar>::value, FieldScalar>
vtmv(const Dune::DenseVector<V1>& v1,
     const FieldScalar m,
     const Dune::DenseVector<V2>& v2)
{
    return m*(v1*v2);
}

} // end namespace Dumux

#endif
