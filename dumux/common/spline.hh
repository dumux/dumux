// $Id$
/*****************************************************************************
 *   Copyright (C) 2008-2010 by Andreas Lauser                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief Provides 3rd order polynomial splines.
 */
#ifndef DUMUX_SPLINE_HH
#define DUMUX_SPLINE_HH

#include "fixedlengthspline_.hh"
#include "variablelengthspline_.hh"
#include "splinecommon_.hh"

namespace Dumux
{
/*!
 * \brief A 3rd order polynomial spline.

 * This class implements a spline \f$s(x)\f$ for which, given \f$n\f$ sampling
 * points \f$x_1, \dots, x_n\f$, the following conditions hold
 *\f{align*}{
   s(x_i)   & = y_i \quad \forall i \in \{1, \dots, n \}\\
   s'(x_1)  & = m_1 \\
   s'(x_n)  & = m_n
   \f}
*
* for any given boundary slopes \f$m_1\f$ and \f$m_n\f$. Alternatively, natural
* splines are supported which are defined by
*\f{align*}{
    s(x_i)     & = y_i \quad \forall i \in \{1, \dots, n \} \\
    s''(x_1)   & = 0 \\
    s''(x_n)   & = 0
\f}
 */


template<class Scalar, int numSamples = 2>
class Spline : public FixedLengthSpline_<Scalar, numSamples>
{
public:
    /*!
     * \brief Default constructor for a spline.
     *
     * To specfiy the acutal curve, use one of the set() methods.
     */
    Spline()
    { };

    /*!
     * \brief Convenience constructor for a natural spline
     *
     * \param x An array containing the \f$x\f$ values of the spline's sampling points
     * \param y An array containing the \f$y\f$ values of the spline's sampling points
     */
    template <class ScalarContainer>
    Spline(const ScalarContainer &x,
           const ScalarContainer &y)
    { this->set(x, y); }

    /*!
     * \brief Convenience constructor for a natural spline
     *
     * \param tuples An array of \f$(x,y)\f$ tuples of the spline's sampling points
     */
    template <class XYContainer>
    Spline(const XYContainer &tuples)
    { this->set(tuples); }

    /*!
     * \brief Convenience constructor for a full spline
     *
     * \param x An array containing the \f$x\f$ values of the spline's sampling points
     * \param y An array containing the \f$y\f$ values of the spline's sampling points
     * \param m0 The slope of the spline at \f$x_0\f$
     * \param m1 The slope of the spline at \f$x_n\f$
     */
    template <class ScalarContainer>
    Spline(const ScalarContainer &x,
           const ScalarContainer &y,
           Scalar m0,
           Scalar m1)
    { this->set(x, y, m0, m1); }

    /*!
     * \brief Convenience constructor for a full spline
     *
     * \param points An array of \f$(x,y)\f$ tuples of the spline's sampling points
     * \param m0 The slope of the spline at \f$x_0\f$
     * \param m1 The slope of the spline at \f$x_n\f$
     */
    template <class XYContainer>
    Spline(const XYContainer &points,
           Scalar m0,
           Scalar m1)
    { this->set(points, m0, m1); }
};

/*!
 * \brief Specialization of a spline with the number of sampling
 *        points only known at run time.
 *
 * This class implements a spline \f$s(x)\f$ for which, given \f$n\f$ sampling
 * points \f$x_1, \dots, x_n\f$, the following conditions hold
 *\f{align*}{
     s(x_i)   & = y_i \quad \forall i \in \{1, \dots, n \}\\
     s'(x_1)  & = m_1 \\
     s'(x_n)  & = m_n
   \f}
 *
 * for any given boundary slopes \f$m_1\f$ and \f$m_n\f$. Alternatively, natural
 * splines are supported which are defined by
 *\f{align*}{
    s(x_i)     & = y_i \quad \forall i \in \{1, \dots, n \} \\
    s''(x_1)   & = 0 \\
    s''(x_n)   & = 0
 \f}
*/
template<class Scalar>
class Spline<Scalar, -1> : public VariableLengthSpline_<Scalar>
{
public:
    /*!
     * \brief Default constructor for a spline.
     *
     * To specfiy the acutal curve, use one of the set() methods.
     */
    Spline()
    { }

    /*!
     * \brief Convenience constructor for a natural spline
     *
     * \param nSamples The number of sampling points (must be > 2)
     * \param x An array containing the \f$x\f$ values of the spline's sampling points
     * \param y An array containing the \f$y\f$ values of the spline's sampling points
     */
    template <class ScalarContainer>
    Spline(int nSamples,
           const ScalarContainer &x,
           const ScalarContainer &y)
    { this->set(nSamples, x, y); }

    /*!
     * \brief Convenience constructor for a natural spline
     *
     * \param nSamples The number of sampling points (must be > 2)
     * \param tuples An array of \f$(x,y)\f$ tuples of the spline's sampling points
     */
    template <class XYContainer>
    Spline(int nSamples,
           const XYContainer &tuples)
    { this->set(nSamples, tuples); }

    /*!
     * \brief Convenience constructor for a natural spline
     *
     * \param x An array containing the \f$x\f$ values of the spline's sampling points (must have a size() method)
     * \param y An array containing the \f$y\f$ values of the spline's sampling points (must have a size() method)
     */
    template <class ScalarContainer>
    Spline(const ScalarContainer &x,
           const ScalarContainer &y)
    { this->set(x, y); }

    /*!
     * \brief Convenience constructor for a natural spline
     *
     * \param tuples An array of \f$(x,y)\f$ tuples of the spline's sampling points (must have a size() method)
     */
    template <class XYContainer>
    Spline(const XYContainer &tuples)
    { this->set(tuples); }

    /*!
     * \brief Convenience constructor for a full spline
     *
     * \param nSamples The number of sampling points (must be >= 2)
     * \param x An array containing the \f$x\f$ values of the spline's sampling points
     * \param y An array containing the \f$y\f$ values of the spline's sampling points
     * \param m0 The slope of the spline at \f$x_0\f$
     * \param m1 The slope of the spline at \f$x_n\f$
     */
    template <class ScalarContainer>
    Spline(int nSamples,
           const ScalarContainer &x,
           const ScalarContainer &y,
           Scalar m0,
           Scalar m1)
    { this->set(nSamples, x, y, m0, m1); }

    /*!
     * \brief Convenience constructor for a full spline
     *
     * \param nSamples The number of sampling points (must be >= 2)
     * \param points An array of \f$(x,y)\f$ tuples of the spline's sampling points
     * \param m0 The slope of the spline at \f$x_0\f$
     * \param m1 The slope of the spline at \f$x_n\f$
     */
    template <class XYContainer>
    Spline(int nSamples,
           const XYContainer &points,
           Scalar m0,
           Scalar m1)
    { this->set(nSamples, points, m0, m1); }

    /*!
     * \brief Convenience constructor for a full spline
     *
     * \param x An array containing the \f$x\f$ values of the spline's sampling points (must have a size() method)
     * \param y An array containing the \f$y\f$ values of the spline's sampling points (must have a size() method)
     * \param m0 The slope of the spline at \f$x_0\f$
     * \param m1 The slope of the spline at \f$x_n\f$
     */
    template <class ScalarContainer>
    Spline(const ScalarContainer &x,
           const ScalarContainer &y,
           Scalar m0,
           Scalar m1)
    { this->set(x, y, m0, m1); }

    /*!
     * \brief Convenience constructor for a full spline
     *
     * \param points An array of \f$(x,y)\f$ tuples of the spline's sampling points (must have a size() method)
     * \param m0 The slope of the spline at \f$x_0\f$
     * \param m1 The slope of the spline at \f$x_n\f$
     */
    template <class XYContainer>
    Spline(const XYContainer &points,
           Scalar m0,
           Scalar m1)
    { this->set(points, m0, m1); }
};

/*!
 * \brief Do not allow splines with zero sampling points
 */
template<class Scalar>
class Spline<Scalar, 0>
// Splines with zero sampling points do not make sense!
{ private: Spline() { }; };

/*!
 * \brief Do not allow splines with one sampling point
 */
template<class Scalar>
class Spline<Scalar, 1>
// Splines with one sampling point do not make sense!
{ private: Spline() { }; };

/*!
 * \brief Spline for two sampling points.
 *
 * For this type of spline there is no natural spline.
 */
template<class Scalar>
class Spline<Scalar, 2> : public SplineCommon_<Scalar, Spline<Scalar, 2> >
{
    friend class  SplineCommon_<Scalar, Spline<Scalar, 2> >;
    typedef Dune::FieldVector<Scalar, 2> Vector;
    typedef Dune::FieldMatrix<Scalar, 2, 2> Matrix;

public:
    Spline()
    {};

    /*!
     * \brief Convenience constructor for a full spline
     *
     * \param x An array containing the \f$x\f$ values of the spline's sampling points
     * \param y An array containing the \f$y\f$ values of the spline's sampling points
     * \param m0 The slope of the spline at \f$x_0\f$
     * \param m1 The slope of the spline at \f$x_n\f$
     */
    template <class ScalarContainer>
    Spline(const ScalarContainer &x,
           const ScalarContainer &y,
           Scalar m0, Scalar m1)
    {
        set(x,y,m0,m1);
    }

    /*!
     * \brief Convenience constructor for a full spline
     *
     * \param points An array of \f$(x,y)\f$ tuples of the spline's sampling points
     * \param m0 The slope of the spline at \f$x_0\f$
     * \param m1 The slope of the spline at \f$x_n\f$
     */
    template <class XYContainer>
    Spline(const XYContainer &points,
           Scalar m0,
           Scalar m1)
    { this->set(points, m0, m1); }

    /*!
     * \brief Convenience constructor for a full spline
     *
     * \param x0 The \f$x\f$ value of the first sampling point
     * \param x1 The \f$x\f$ value of the second sampling point
     * \param y0 The \f$y\f$ value of the first sampling point
     * \param y1 The \f$y\f$ value of the second sampling point
     * \param m0 The slope of the spline at \f$x_0\f$
     * \param m1 The slope of the spline at \f$x_n\f$
     */
    Spline(Scalar x0, Scalar x1,
           Scalar y0, Scalar y1,
           Scalar m0, Scalar m1)
    {
        set(x0, x1,
            y0, y1,
            m0, m1);
    };

    /*!
     * \brief Returns the number of sampling points.
     */
    int numSamples() const
    { return 2; }

    /*!
     * \brief Set the sampling points and the boundary slopes of the
     *        spline.
     *
     * \param x0 The \f$x\f$ value of the first sampling point
     * \param x1 The \f$x\f$ value of the second sampling point
     * \param y0 The \f$y\f$ value of the first sampling point
     * \param y1 The \f$y\f$ value of the second sampling point
     * \param m0 The slope of the spline at \f$x_0\f$
     * \param m1 The slope of the spline at \f$x_1\f$
     */
    void set(Scalar x0, Scalar x1,
             Scalar y0, Scalar y1,
             Scalar m0, Scalar m1)
    {
        Matrix M(numSamples());
        Vector d;
        assignXY_(x0, x1, y0, y1);
        this->makeFullSystem_(M, d, m0, m1);

        // solve for the moments
        M.solve(m_, d);
    }

    /*!
     * \brief Set the sampling points and the boundary slopes of the
     *        spline.
     *
     * \param x An array containing the \f$x\f$ values of the sampling points
     * \param y An array containing the \f$y\f$ values of the sampling points
     * \param m0 The slope of the spline at \f$x_0\f$
     * \param m1 The slope of the spline at \f$x_1\f$
     */
    template <class ScalarContainer>
    void set(const ScalarContainer &x, const ScalarContainer &y,
             Scalar m0, Scalar m1)
    {
        Matrix M(numSamples());
        Vector d;
        assignXY_(x[0], x[1], y[0], y[1]);
        this->makeFullSystem_(M, d, m0, m1);

        // solve for the moments
        M.solve(m_, d);
    }

    /*!
     * \brief Set the sampling points and the boundary slopes of the
     *        spline.
     *
     * \param points An array of \f$(x,y)\f$ tuples of the spline's sampling points
     * \param m0 The slope of the spline at \f$x_0\f$
     * \param m1 The slope of the spline at \f$x_1\f$
     */
    template <class XYContainer>
    void set(const XYContainer &points, Scalar m0, Scalar m1)
    {
        Matrix M;
        Vector d;
        assignXY_(points[0][0], points[1][0], points[0][1], points[1][1]);
        this->makeFullSystem_(M, d, m0, m1);

        // solve for the moments
        M.solve(m_, d);
    }

protected:
    void assignXY_(Scalar x0, Scalar x1,
                   Scalar y0, Scalar y1)
    {
        if (x0 > x1) {
            xPos_[0] = x1;
            xPos_[1] = x0;
            yPos_[0] = y1;
            yPos_[1] = y0;
        }
        else {
            xPos_[0] = x0;
            xPos_[1] = x1;
            yPos_[0] = y0;
            yPos_[1] = y1;
        }
    };

    /*!
     * \brief Returns the x coordinate of the i-th sampling point.
     */
    Scalar x_(int i) const
    { return xPos_[i]; }

    /*!
     * \brief Returns the y coordinate of the i-th sampling point.
     */
    Scalar y_(int i) const
    { return yPos_[i]; }

    /*!
     * \brief Returns the moment (i.e. second derivative) of the
     *        spline at the i-th sampling point.
     */
    Scalar moment_(int i) const
    { return m_[i]; }

    Vector xPos_;
    Vector yPos_;
    Vector m_;
};

}

#endif
