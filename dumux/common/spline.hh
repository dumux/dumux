// $Id$
/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 * \brief A 3rd order polynomial p(x) for which, given two points x1 and x2!=x1,
 *        the following hold:
 *        p(x1) = y1
 *        p(x2) = y2
 *        p'(x1) = m1
 *        p'(x2) = m2
 *
 * or any given y1, y2, m1, m2.
 */
#ifndef DUMUX_SPLINE_HH
#define DUMUX_SPLINE_HH

#include <dumux/common/valgrind.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/istl/btdmatrix.hh>

namespace Dumux
{

/*
 * \brief A 3rd order polynomial spline.

 * This class implements a spline s(x) for which, given n sampling
 * points x_1, ..., x_n, the following holds:
 *
 *        s(x_i) = y_i
 *        s'(x_1) = m_1
 *        s'(x_n) = m_n
 *
 * for any given boundary slopes m_1 and m_n.
 */
template<class ScalarT, int numSamples = 2>
class Spline
{
    typedef ScalarT Scalar;
    typedef Dune::FieldVector<Scalar,numSamples> FieldVector;

public:
    Spline()
    { Valgrind::SetUndefined(*this); };

    Spline(const FieldVector &x,
           const FieldVector &y,
           Scalar m0,
           Scalar m1)
    {
        set(x, y, m0, m1);
    }

    Spline(const Scalar *x,
           const Scalar *y,
           Scalar m0,
           Scalar m1)
    {
        set(x, y, m0, m1);
    }

    /*!
     * \brief Set the sampling points and the boundary slopes of the
     *        spline function.
     */
    void set(const FieldVector &x,
             const FieldVector &y,
             Scalar m0,
             Scalar m1)
    {
        x_ = x;
        y_ = y;

        // See: J. Stoer: "Numerische Mathematik 1", 9th edition,
        // Springer, 2005, p. 111
        Dune::FieldMatrix<Scalar, numSamples, numSamples> M(0);
        FieldVector d(0);

        const int n = numSamples - 1;

        // second to next to last rows
        for (int i = 1; i < n; ++i) {
            const Scalar lambda_i = h_(i + 1) / (h_(i) + h_(i+1));
            const Scalar mu_i = 1 - lambda_i;
            const Scalar d_i =
                6 / (h_(i) + h_(i+1))
                *
                ( (y_[i+1] - y_[i])/h_(i+1) - (y_[i] - y_[i-1])/h_(i));

            M[i][i-1] = mu_i;
            M[i][i] = 2;
            M[i][i + 1] = lambda_i;
            d[i] = d_i;
        };

        // first row
        M[0][0] = 2;
        M[0][1] = 1;
        d[0] = 6/h_(1) * ( (y_[1] - y_[0])/h_(1) - m0);

        // last row
        M[n][n] = 2;
        M[n][n - 1] = 1;
        d[n] =
            6/h_(n)
            *
            (m1 - (y_[n] - y_[n - 1])/h_(n));

        // solve for the moments
        M.solve(moments_, d);
    }

    /*!
     * \brief Set the sampling points and the boundary slopes of the
     *        spline function.
     */
    void set(const Scalar *x,
             const Scalar *y,
             Scalar m0,
             Scalar m1)
    {
        set(toFieldVector_(x),
            toFieldVector_(y),
            m0, m1);
    };

    /*!
     * \brief Return true iff the given x is in range [x1, xn].
     */
    bool applies(Scalar x) const
    { return x_[0] <= x && x <= x_[numSamples - 1]; };

    /*!
     * \brief Return the x value of the leftmost sampling point.
     */
    Scalar xMin() const
    { return x_[0]; };

    /*!
     * \brief Return the x value of the rightmost sampling point.
     */
    Scalar xMax() const
    { return x_[numSamples - 1]; };

    /*!
     * \brief Evaluate the spline at a given position.
     */
    Scalar eval(Scalar x) const
    {
        assert(applies(x));

        // See: J. Stoer: "Numerische Mathematik 1", 9th edition,
        // Springer, 2005, p. 109

        int i = findIntervalIdx_(x);

        Scalar h_i1 = h_(i + 1);
        Scalar x_i = x - x_[i];
        Scalar x_i1 = x_[i+1] - x;

        Scalar A_i =
            (y_[i+1] - y_[i])/h_i1
            -
            h_i1/6*(moments_[i+1] - moments_[i]);
        Scalar B_i = y_[i] - moments_[i]* (h_i1*h_i1) / 6;

        return
            moments_[i]* x_i1*x_i1*x_i1 / (6 * h_i1)
            +
            moments_[i + 1]* x_i*x_i*x_i / (6 * h_i1)
            +
            A_i*x_i
            +
            B_i;
    }

    /*!
     * \brief Evaluate the spline's derivative at a given position.
     */
    Scalar evalDerivative(Scalar x) const
    {
        assert(applies(x));

        // See: J. Stoer: "Numerische Mathematik 1", 9th edition,
        // Springer, 2005, p. 109

        int i = findIntervalIdx_(x);

        Scalar h_i1 = h_(i + 1);
        Scalar x_i = x - x_[i];
        Scalar x_i1 = x_[i+1] - x;

        Scalar A_i =
            (y_[i+1] - y_[i])/h_i1
            -
            h_i1/6*(moments_[i+1] - moments_[i]);

        return
            -moments_[i] * x_i1*x_i1 / (2 * h_i1)
            +
            moments_[i + 1] * x_i*x_i / (2 * h_i1)
            +
            A_i;
    }

    /*!
     * \brief Returns 1 if the spline is monotonically increasing, -1
     *        if the spline is mononously decreasing and 0 if the
     *        spline is not monotonous in the interval (x0, x1).
     *
     * In the corner case where the whole spline is flat, it returns
     * 2.
     */
    int monotonic(Scalar x0, Scalar x1) const
    {
        assert(applies(x0));
        assert(applies(x1));
        assert(x0 <= x1);

        // corner case where the whole spline is a constant
        if (moments_[0] == 0 &&
            moments_[1] == 0 &&
            y_[0] == y_[1])
        {
            // actually the is monotonically increasing as well as
            // monotonously decreasing
            return 2;
        }


        int i = findIntervalIdx_(x0);
        if (x_[i] <= x0 && x1 <= x_[i+1]) {
            // the interval to check is completely included in a
            // single segment
            return monotonic_(i, x0, x1);
        }

        // make sure that the segments which are completly in the
        // interval [x0, x1] all exhibit the same monotonicity.
        int iEnd = findIntervalIdx_(x1);
        int r = monotonic_(i, x0, x_[1]);
        for (; i < iEnd - 1; ++i)
            if (r != monotonic_(i, x_[i], x_[i + 1]))
                return 0;

        // check for the last segment
        if (x_[iEnd] < x1 && r != monotonic_(iEnd, x_[iEnd], x1))
        { return 0; }

        return r;
    }

    /*!
     * \brief Same as monotonic(x0, x1), but with the entire range of the
     *        spline as interval.
     */
    int monotonic() const
    { return monotonic(x_[0], x_[numSamples - 1]); }

    /*!
     * \brief Prints k tuples of the format (x, y, dx/dy, isMonotonic)
     *        to stdout.
     *
     * If the spline does not apply for parts of [x0, x1] it is
     * extrapolated using a straight line. The result can be inspected
     * using the following commands:
     *
     ----------- snip -----------
     ./yourProgramm > spline.csv
     gnuplot

     gnuplot> plot "spline.csv" using 1:2 w l ti "Curve", \
     "spline.csv" using 1:3 w l ti "Derivative", \
     "spline.csv" using 1:4 w p ti "Monotonic"
     ----------- snap -----------
     */
    void printCSV(Scalar x0, Scalar x1, int k) const
    {
        const int n = numSamples - 1;
        for (int i = 0; i <= k; ++i) {
            double x = i*(x1 - x0)/k + x0;
            double x_p1 = x + (x1 - x0)/k;
            double y;
            double dy_dx;
            double mono = 1;
            if (!applies(x)) {
                if (x < 0) {
                    dy_dx = evalDerivative(x_[0]);
                    y = (x - x_[0])*dy_dx + y_[0];
                    mono = (dy_dx>0)?1:-1;
                }
                else if (x > x_[n]) {
                    dy_dx = evalDerivative(x_[n]);
                    y = (x - x_[n])*dy_dx + y_[n];
                    mono = (dy_dx>0)?1:-1;
                }
                else {
                    std::cerr << "ooops: " << x << "\n";
                    exit(1);
                }
            }
            else {
                y = eval(x);
                dy_dx = evalDerivative(x);
                mono = monotonic(std::max(x_[0], x), std::min(x_[n], x_p1));
            }

            std::cout << x << " " << y << " " << dy_dx << " " << mono << "\n";
        }
    }

private:
    // returns the monotonicality of an interval of a spline segment
    int monotonic_(int i, Scalar x0, Scalar x1) const
    {
        // shift the interval so that it is consistent with the
        // definitions by Stoer
        x0 = x0 - x_[i];
        x1 = x1 - x_[i];

        Scalar a = a_(i);
        Scalar b = b_(i);
        Scalar c = c_(i);

        Scalar disc = 4*b*b - 12*a*c;
        if (disc < 0) {
            // discriminant is smaller than 0, i.e. the segment does
            // not exhibit any extrema.
            return (x0*(x0*3*a + 2*b) + c > 0) ? 1 : -1;
        }
        disc = std::sqrt(disc);
        Scalar xE1 = (-2*b + disc)/(6*a);
        Scalar xE2 = (-2*b - disc)/(6*a);

        if (disc == 0) {
            // saddle point -> no extrema
            if (xE1 == x0)
                // make sure that we're not picking the saddle point
                // to determine whether we're monotonically increasing
                // or decreasing
                x0 = x1;
            return (x0*(x0*3*a + 2*b) + c > 0) ? 1 : -1;
        };
        if ((x0 < xE1 && xE1 < x1) ||
            (x0 < xE2 && xE2 < x1))
        {
            // there is an extremum in the range (x0, x1)
            return 0;
        }
        // no extremum in range (x0, x1)
        x0 = (x0 + x1)/2; // pick point in the middle of the interval
                          // to avoid extrema on the boundaries
        return (x0*(x0*3*a + 2*b) + c > 0) ? 1 : -1;
    };

    Scalar h_(int i) const
    { return x_[i] - x_[i-1]; }

    int findIntervalIdx_(Scalar x) const
    {
        // bisection
        int iLow = 0;
        int iHigh = numSamples - 1;

        while (iLow + 1 < iHigh) {
            int i = (iLow + iHigh) / 2;
            if (x_[i] > x)
                iHigh = i;
            else
                iLow = i;
        };
        return iLow;
    };

    // returns the coefficient in front of the x^3 term. In Stoer this
    // is delta.
    Scalar a_(int i) const
    { return (moments_[i+1] - moments_[i])/(6*h_(i+1)); }

    // returns the coefficient in front of the x^2 term In Stoer this
    // is gamma.
    Scalar b_(int i) const
    { return moments_[i]/2; }

    // returns the coefficient in front of the x term. In Stoer this
    // is beta.
    Scalar c_(int i) const
    { return (y_[i+1] - y_[i])/h_(i + 1) - h_(i+1)/6*(2*moments_[i] + moments_[i+1]); }

    // returns the coefficient in front of the constant term. In Stoer this
    // is alpha.
    Scalar d_(int i) const
    { return y_[i]; }

    const FieldVector &toFieldVector_(const Scalar *v) const
    { return *reinterpret_cast<const FieldVector*>(v); };

    FieldVector moments_; // moments
    FieldVector x_; // x positions of the sampling points
    FieldVector y_; // y positions of the sampling points
};

/*
 * \brief A 3rd order polynomial spline.

 * This class implements a spline s(x) for which, given n (known at
 * run-time) sampling points x_1, ..., x_n, the following holds:
 *
 *        s(x_i) = y_i
 *        s'(x_1) = m_1
 *        s'(x_n) = m_n
 *
 * for any given boundary slopes m_1 and m_n.
 */
template<class ScalarT>
class Spline<ScalarT, -1>
{
    typedef ScalarT Scalar;
    typedef Dune::BlockVector<Dune::FieldVector<Scalar,1> > BlockVector;
    typedef Dune::BTDMatrix<Dune::FieldMatrix<Scalar, 1, 1> > BTDMatrix;

public:
    Spline()
    { };

    /*!
     * \brief Set the sampling points and the boundary slopes of the
     *        spline function.
     */
    void set(int numSamples,
             const Scalar *x,
             const Scalar *y,
             Scalar m0,
             Scalar m1)
    {
        BTDMatrix M(numSamples);
        BlockVector d(numSamples);
        fillNatural_(M, d, numSamples, x, y);

        int n = numSamples - 1;
        // first row
        M[0][1] = 1;
        d[0] = 6/h_(1) * ( (y_[1] - y_[0])/h_(1) - m0);

        // last row
        M[n][n - 1] = 1;
        d[n] =
            6/h_(n)
            *
            (m1 - (y_[n] - y_[n - 1])/h_(n));

        // solve for the moments
        M.solve(moments_, d);
    }

    /*!
     * \brief Set the sampling points  of the
     *        spline function.
     *
     * The second derivatives at the boundaries are 0 in this
     * case. (i.e. this is a natural spline)
     */
    void set(int numSamples,
             const Scalar *x,
             const Scalar *y)
    {
        BTDMatrix M(numSamples);
        BlockVector d(numSamples);
        fillNatural_(M, d, numSamples, x, y);

        // solve for the moments
        M.solve(moments_, d);
    }

    /*!
     * \brief Set the sampling points  of the
     *        spline function.
     *
     * The second derivatives at the boundaries are 0 in this
     * case. (i.e. this is a natural spline)
     */
    void set(int numSamples,
             const BlockVector &x,
             const BlockVector &y)
    {
        BTDMatrix M(numSamples);
        BlockVector d(numSamples);
        fillNatural_(M, d, numSamples, x, y);

        // solve for the moments
        M.solve(moments_, d);
    }

    /*!
     * \brief Return true iff the given x is in range [x1, xn].
     */
    bool applies(Scalar x) const
    {
        return x_[0] <= x && x <= x_[numSamples() - 1];
    };

    /*!
     * \brief Return the x value of the leftmost sampling point.
     */
    Scalar xMin() const
    { return x_[0]; };

    /*!
     * \brief Return the x value of the rightmost sampling point.
     */
    Scalar xMax() const
    { return x_[numSamples() - 1]; };

    /*!
     * \brief Returns the number of sampling points.
     */
    int numSamples() const
    { return x_.size(); }

    /*!
     * \brief Evaluate the spline at a given position.
     */
    Scalar eval(Scalar x) const
    {
        assert(applies(x));

        // See: J. Stoer: "Numerische Mathematik 1", 9th edition,
        // Springer, 2005, p. 109

        int i = findIntervalIdx_(x);

        Scalar h_i1 = h_(i + 1);
        Scalar x_i = x - x_[i];
        Scalar x_i1 = x_[i+1] - x;

        Scalar A_i =
            (y_[i+1] - y_[i])/h_i1
            -
            h_i1/6*(moments_[i+1] - moments_[i]);
        Scalar B_i = y_[i] - moments_[i]* (h_i1*h_i1) / 6;

        return
            moments_[i]* x_i1*x_i1*x_i1 / (6 * h_i1)
            +
            moments_[i + 1]* x_i*x_i*x_i / (6 * h_i1)
            +
            A_i*x_i
            +
            B_i;
    }

    /*!
     * \brief Evaluate the spline's derivative at a given position.
     */
    Scalar evalDerivative(Scalar x) const
    {
        assert(applies(x));

        // See: J. Stoer: "Numerische Mathematik 1", 9th edition,
        // Springer, 2005, p. 109

        int i = findIntervalIdx_(x);

        Scalar h_i1 = h_(i + 1);
        Scalar x_i = x - x_[i];
        Scalar x_i1 = x_[i+1] - x;

        Scalar A_i =
            (y_[i+1] - y_[i])/h_i1
            -
            h_i1/6*(moments_[i+1] - moments_[i]);

        return
            -moments_[i] * x_i1*x_i1 / (2 * h_i1)
            +
            moments_[i + 1] * x_i*x_i / (2 * h_i1)
            +
            A_i;
    }

    /*!
     * \brief Returns 1 if the spline is monotonically increasing, -1
     *        if the spline is mononously decreasing and 0 if the
     *        spline is not monotonous in the interval (x0, x1).
     *
     * In the corner case where the whole spline is flat, it returns
     * 2.
     */
    int monotonic(Scalar xi0, Scalar xi1) const
    {
        assert(applies(xi0));
        assert(applies(xi1));

        Scalar x0 = std::min(xi0, xi1);
        Scalar x1 = std::max(xi0, xi1);

        // corner case where the whole spline is a constant
        if (moments_[0] == 0 &&
            moments_[1] == 0 &&
            y_[0] == y_[1])
        {
            // actually the is monotonically increasing as well as
            // monotonously decreasing
            return 2;
        }


        int i = findIntervalIdx_(x0);
        if (x_[i] <= x0 && x1 <= x_[i+1]) {
            // the interval to check is completely included in a
            // single segment
            return monotonic_(i, x0, x1);
        }

        // make sure that the segments which are completly in the
        // interval [x0, x1] all exhibit the same monotonicity.
        int iEnd = findIntervalIdx_(x1);
        int r = monotonic_(i, x0, x_[1]);
        for (; i < iEnd - 1; ++i)
            if (r != monotonic_(i, x_[i], x_[i + 1]))
                return 0;

        // check for the last segment
        if (x_[iEnd] < x1 && r != monotonic_(iEnd, x_[iEnd], x1))
        { return 0; }

        return r;
    }

    /*!
     * \brief Same as monotonic(x0, x1), but with the entire range of the
     *        spline as interval.
     */
    int monotonic() const
    { return monotonic(x_[0], x_[numSamples() - 1]); }

    /*!
     * \brief Prints k tuples of the format (x, y, dx/dy, isMonotonic)
     *        to stdout.
     *
     * If the spline does not apply for parts of [x0, x1] it is
     * extrapolated using a straight line. The result can be inspected
     * using the following commands:
     *
     ----------- snip -----------
     ./yourProgramm > spline.csv
     gnuplot

     gnuplot> plot "spline.csv" using 1:2 w l ti "Curve", \
     "spline.csv" using 1:3 w l ti "Derivative", \
     "spline.csv" using 1:4 w p ti "Monotonic"
     ----------- snap -----------
     */
    void printCSV(Scalar xi0, Scalar xi1, int k) const
    {
        Scalar x0 = std::min(xi0, xi1);
        Scalar x1 = std::max(xi0, xi1);
        const int n = numSamples() - 1;
        for (int i = 0; i <= k; ++i) {
            double x = i*(x1 - x0)/k + x0;
            double x_p1 = x + (x1 - x0)/k;
            double y;
            double dy_dx;
            double mono = 1;
            if (!applies(x)) {
                if (x < 0) {
                    dy_dx = evalDerivative(x_[0]);
                    y = (x - x_[0])*dy_dx + y_[0];
                    mono = (dy_dx>0)?1:-1;
                }
                else if (x > x_[n]) {
                    dy_dx = evalDerivative(x_[n]);
                    y = (x - x_[n])*dy_dx + y_[n];
                    mono = (dy_dx>0)?1:-1;
                }
                else {
                    std::cerr << "ooops: " << x << "\n";
                    exit(1);
                }
            }
            else {
                y = eval(x);
                dy_dx = evalDerivative(x);
                mono = monotonic(std::max<Scalar>(x_[0][0], x), std::min<Scalar>(x_[n][0], x_p1));
            }

            std::cout << x << " " << y << " " << dy_dx << " " << mono << "\n";
        }
    }

private:
    // fill the system of equations for a natrural spline
    template <class ScalarArray>
    void fillNatural_(BTDMatrix &M,
                      BlockVector &d,
                      int numSamples,
                      const ScalarArray &x,
                      const ScalarArray &y)
    {
        x_.resize(numSamples);
        y_.resize(numSamples);
        moments_.resize(numSamples);

        // copy sample points, make sure that the first x value is
        // smaller than the last one
        for (int i = 0; i < numSamples; ++i) {
            int idx = i;
            if (x[0] > x[numSamples - 1])
                idx = numSamples - i - 1;
            x_[i] = x[idx];
            y_[i] = y[idx];
        }

        M = 0;
        d = 0;

        // See: J. Stoer: "Numerische Mathematik 1", 9th edition,
        // Springer, 2005, p. 111
        const int n = numSamples - 1;

        // second to next to last rows
        for (int i = 1; i < n; ++i) {
            const Scalar lambda_i = h_(i + 1) / (h_(i) + h_(i+1));
            const Scalar mu_i = 1 - lambda_i;
            const Scalar d_i =
                6 / (h_(i) + h_(i+1))
                *
                ( (y_[i+1] - y_[i])/h_(i+1) - (y_[i] - y_[i-1])/h_(i));

            M[i][i-1] = mu_i;
            M[i][i] = 2;
            M[i][i + 1] = lambda_i;
            d[i] = d_i;
        };

        // first row
        M[0][0] = 2;

        // last row
        M[n][n] = 2;
    }

    // returns the monotonicality of an interval of a spline segment
    int monotonic_(int i, Scalar x0, Scalar x1) const
    {
        // shift the interval so that it is consistent with the
        // definitions by Stoer
        x0 = x0 - x_[i];
        x1 = x1 - x_[i];

        Scalar a = a_(i);
        Scalar b = b_(i);
        Scalar c = c_(i);

        Scalar disc = 4*b*b - 12*a*c;
        if (disc < 0) {
            // discriminant is smaller than 0, i.e. the segment does
            // not exhibit any extrema.
            return (x0*(x0*3*a + 2*b) + c > 0) ? 1 : -1;
        }
        disc = std::sqrt(disc);
        Scalar xE1 = (-2*b + disc)/(6*a);
        Scalar xE2 = (-2*b - disc)/(6*a);

        if (disc == 0) {
            // saddle point -> no extrema
            if (xE1 == x0)
                // make sure that we're not picking the saddle point
                // to determine whether we're monotonically increasing
                // or decreasing
                x0 = x1;
            return (x0*(x0*3*a + 2*b) + c > 0) ? 1 : -1;
        };
        if ((x0 < xE1 && xE1 < x1) ||
            (x0 < xE2 && xE2 < x1))
        {
            // there is an extremum in the range (x0, x1)
            return 0;
        }
        // no extremum in range (x0, x1)
        x0 = (x0 + x1)/2; // pick point in the middle of the interval
                          // to avoid extrema on the boundaries
        return (x0*(x0*3*a + 2*b) + c > 0) ? 1 : -1;
    };

    Scalar h_(int i) const
    {
        assert(x_[i] - x_[i-1] > 0);

        return x_[i] - x_[i-1];
    }

    int findIntervalIdx_(Scalar x) const
    {
        // bisection
        int iLow = 0;
        int iHigh = numSamples() - 1;

        while (iLow + 1 < iHigh) {
            int i = (iLow + iHigh) / 2;
            if (x_[i] > x)
                iHigh = i;
            else
                iLow = i;
        };
        return iLow;
    };

    // returns the coefficient in front of the x^3 term. In Stoer this
    // is delta.
    Scalar a_(int i) const
    { return (moments_[i+1] - moments_[i])/(6*h_(i+1)); }

    // returns the coefficient in front of the x^2 term In Stoer this
    // is gamma.
    Scalar b_(int i) const
    { return moments_[i]/2; }

    // returns the coefficient in front of the x term. In Stoer this
    // is beta.
    Scalar c_(int i) const
    { return (y_[i+1] - y_[i])/h_(i + 1) - h_(i+1)/6*(2*moments_[i] + moments_[i+1]); }

    // returns the coefficient in front of the constant term. In Stoer this
    // is alpha.
    Scalar d_(int i) const
    { return y_[i]; }

    const BlockVector &toBlockVector_(const Scalar *v) const
    { return *reinterpret_cast<const BlockVector*>(v); };

    BlockVector moments_; // moments
    BlockVector x_; // x positions of the sampling points
    BlockVector y_; // y positions of the sampling points
};

/*
 * \brief A 3rd order polynomial p(x) = a x^3 + b x^2 + c x + d for
 *        which, given two distinct points x1 and x2, the following holds:
 *        p(x1) = y1
 *        p(x2) = y2
 *        p'(x1) = m1
 *        p'(x2) = m2
 *
 * for any given y1, y2, m1, m2.
 */
template<class ScalarT>
class Spline<ScalarT, 2>
{
    typedef ScalarT Scalar;
    typedef Dune::FieldVector<Scalar,2> FieldVector;

public:
    Spline()
    {
        Valgrind::SetUndefined(*this);
    };

    Spline(Scalar x1,
           Scalar x2,
           Scalar y1,
           Scalar y2,
           Scalar m1,
           Scalar m2)
    {
        set(x1, x2, y1, y2, m1, m2);
    }

    Spline(const FieldVector &x,
           const FieldVector &y,
           Scalar m1,
           Scalar m2)
    { set(x, y, m1, m2); }

    Spline(const Scalar *x,
           const Scalar *y,
           Scalar m0,
           Scalar m1)
    {
        set(x, y, m0, m1);
    }

    /*!
     * \brief Set the parameters of the spline.
     *
     * This method is the same as in the generic case of more than 2
     * sampling points.
     */
    void set(const FieldVector &x,
             const FieldVector &y,
             Scalar m1,
             Scalar m2)
    { set(x[0], x[1], y[0], y[1], m1, m2); }

    /*!
     * \brief Set the parameters of the spline.
     *
     * This method is the same as in the generic case of more than 2
     * sampling points.
     */
    void set(const Scalar *x,
             const Scalar *y,
             Scalar m1,
             Scalar m2)
    { set(x[0], x[1], y[0], y[1], m1, m2); };

    /*!
     * \brief Set the parameters of the spline.
     *
     * This is a method provided for convinience.
     */
    void set(Scalar x1,
             Scalar x2,
             Scalar y1,
             Scalar y2,
             Scalar m1,
             Scalar m2)
    {
        assert(x1 != x2);

        x1_ = x1;
        x2_ = x2;

        // calculate the coefficents of the polynomial. this
        // is pretty cumbersome, since we basically solve a
        // 4x4 matrix analytically here.
        Scalar tmpLeft  = 2*(x1*x1*x1 - x2*x2*x2) - 3*(x1 - x2)*(x1*x1 + x2*x2);
        Scalar tmpRight = 2*(y1 - y2) - (x1 - x2)*(m1 + m2);
        a_ = tmpRight/tmpLeft;

        tmpLeft = 2*(x1 - x2);
        tmpRight = m1 - m2 - (3*(x1*x1 - x2*x2)*a_);
        b_ = tmpRight/tmpLeft;

        tmpRight = m1 - (3*x1*x1*a_ + 2*x1*b_);
        tmpLeft = 1;
        c_ = tmpRight/tmpLeft;

        tmpRight = y1 - (x1*x1*x1*a_ + x1*x1*b_ + x1*c_);
        tmpLeft = 1;
        d_ = tmpRight/tmpLeft;

        assert(fabs(eval(x1) - y1) < std::max(1e-11, fabs(1e-8*y1)));
        assert(fabs(eval(x2) - y2) < std::max(1e-11, fabs(1e-8*y2)));
        assert(fabs(evalDerivative(x1) - m1) < std::max(1e-11, fabs(1e-8*m1)));
        assert(fabs(evalDerivative(x2) - m2) < std::max(1e-11, fabs(1e-8*m2)));
    }

    /*!
     * \brief Return true iff the given x is in range [x1, x2].
     */
    bool applies(Scalar x) const
    { return x1_ <= x && x <= x2_; };


    /*!
     * \brief Return the x value of the leftmost sampling point.
     */
    Scalar xMin() const
    { return x1_; };

    /*!
     * \brief Return the x value of the rightmost sampling point.
     */
    Scalar xMax() const
    { return x2_; };

    /*!
     * \brief Evaluate the polynomial at a given position.
     */
    Scalar eval(Scalar x) const
    {
        assert(applies(x));
        return ((a_*x + b_)*x + c_)*x + d_;
    }

    /*!
     * \brief Evaluate the polynomial's derivative at a given position.
     */
    Scalar evalDerivative(Scalar x) const
    {
        assert(applies(x));
        return (3*a_*x + 2*b_)*x + c_;
    }

    /*!
     * \brief Returns 1 if the spline is monotonically increasing, -1
     *        if the spline is mononously decreasing and 0 if the
     *        spline is not monotonous in the interval [x0, x1].
     *
     * In the corner case where the whole spline is flat, it returns
     * 2.
     */
    int monotonic(Scalar x0, Scalar x1) const
    {
        // corner case: flat line
        if (a_ == 0 && b_ == 0 && c_ == 0)
            return 2;

        Scalar disc = 4*b_*b_ - 12*a_*c_;
        if (disc < 0) {
            // discriminant is smaller than 0, i.e. the segment does
            // not exhibit any extrema.
            return (x0*(x0*3*a_ + 2*b_) + c_ > 0) ? 1 : -1;
        }

        disc = std::sqrt(disc);
        Scalar xE1 = (-2*b_ + disc)/(6*a_);
        Scalar xE2 = (-2*b_ - disc)/(6*a_);
        if (disc == 0) {
            // saddle point -> no extrema
            if (xE1 == x0)
                // make sure that we're not picking the saddle point
                // to determine whether we're monotonically increasing
                // or decreasing
                x0 = x1;
            return (x0*(x0*3*a_ + 2*b_) + c_ > 0) ? 1 : -1;
        };
        if ((x0 < xE1 && xE1 < x1) ||
            (x0 < xE2 && xE2 < x1))
        {
            // there is an extremum in the range (x0, x1)
            return 0;
        }

        // no extremum in range (x0, x1)
        x0 = (x0 + x1)/2; // pick point in the middle of the interval
                          // to avoid extrema on the boundaries
        return (x0*(x0*3*a_ + 2*b_) + c_ > 0) ? 1 : -1;
    };

    /*!
     * \brief Returns 1 if the spline is monotonically increasing, -1
     *        if the spline is mononously decreasing and 0 if the
     *        spline is not monotonous over the whole interval where
     *        it applies.
     *
     * In the corner case where the whole spline is flat, it returns
     * 2.
     */
    int monotonic() const
    { return monotonic(x1_, x2_); }

    /*!
     * \brief Prints k tuples of the format (x, y, dx/dy, isMonotonic)
     *        to stdout.
     *
     * If the spline does not apply for parts of [x0, x1] it is
     * extrapolated using a straight line. The result can be inspected
     * using the following commands:
     *
     ----------- snip -----------
     ./yourProgramm > spline.csv
     gnuplot

     gnuplot> plot "spline.csv" using 1:2 w l ti "Curve", \
     "spline.csv" using 1:3 w l ti "Derivative", \
     "spline.csv" using 1:4 w p ti "Monotonic"
     ----------- snap -----------
     */
    void printCSV(Scalar x0, Scalar x1, int k) const
    {
        for (int i = 0; i <= k; ++i) {
            double x = i*(x1 - x0)/k + x0;
            double x_p1 = x + (x1 - x0)/k;
            double y;
            double dy_dx;
            double mono = 1;
            if (!applies(x)) {
                if (x < 0) {
                    dy_dx = evalDerivative(x1_);
                    y = (x - x1_)*dy_dx + eval(x1_);
                    mono = (dy_dx>0)?1:-1;
                }
                else if (x > x2_) {
                    dy_dx = evalDerivative(x2_);
                    y = (x - x2_)*dy_dx + eval(x2_);
                    mono = (dy_dx>0)?1:-1;
                }
                else {
                    std::cerr << "ooops: " << x << "\n";
                    exit(1);
                }
            }
            else {
                y = eval(x);
                dy_dx = evalDerivative(x);
                mono = monotonic(std::max(x1_, x), std::min(x2_, x_p1));
            }

            std::cout << x << " " << y << " " << dy_dx << " " << mono << "\n";
        }
    }

private:
    Scalar a_;
    Scalar b_;
    Scalar c_;
    Scalar d_;

    Scalar x1_;
    Scalar x2_;
};

}

#endif
