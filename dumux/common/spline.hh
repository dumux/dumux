// $Id$
/*****************************************************************************
 *   Copyright (C) 2008-2010 by Andreas Lauser                               *
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

#include "fixedlengthspline_.hh"
#include "variablelengthspline_.hh"
#include "splinecommon_.hh"

namespace Dumux
{
/*!
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
template<class Scalar, int numSamples = 2>
class Spline : public FixedLengthSpline_<Scalar, numSamples>
{
public: 
    Spline()
    { };

    // convenience constructor
    template <class ScalarContainer>
    Spline(const ScalarContainer &x, 
           const ScalarContainer &y)
    { this->set(x, y); }
    
    // convenience constructor
    template <class XYContainer>
    Spline(const XYContainer &tuples)
    { this->set(tuples); }

    // convenience constructor
    template <class ScalarContainer>
    Spline(const ScalarContainer &x,
           const ScalarContainer &y,
           Scalar m0,
           Scalar m1)
    { this->set(x, y, m0, m1); }

    // convenience constructor
    template <class XYContainer>
    Spline(const XYContainer &points,
           Scalar m0,
           Scalar m1)
    { this->set(points, m0, m1); }
};

/*
 * \brief Specialization of a spline with the number of sampling
 *        points only known at run time.
 */
template<class Scalar>
class Spline<Scalar, -1> : public VariableLengthSpline_<Scalar>
{ 
public:
    Spline()
    { }

    // convenience constructor
    template <class ScalarContainer>
    Spline(int nSamples,
           const ScalarContainer &x, 
           const ScalarContainer &y)
    { this->set(nSamples, x, y); }
    
    // convenience constructor
    template <class XYContainer>
    Spline(int nSamples,
           const XYContainer &tuples)
    { this->set(nSamples, tuples); }

    // convenience constructor
    template <class ScalarContainer>
    Spline(const ScalarContainer &x, 
           const ScalarContainer &y)
    { this->set(x, y); }
    
    // convenience constructor
    template <class XYContainer>
    Spline(const XYContainer &tuples)
    { this->set(tuples); }

    // convenience constructor
    template <class ScalarContainer>
    Spline(int nSamples, 
           const ScalarContainer &x,
           const ScalarContainer &y,
           Scalar m0,
           Scalar m1)
    { this->set(nSamples, x, y, m0, m1); }

    // convenience constructor
    template <class XYContainer>
    Spline(int nSamples, 
           const XYContainer &points,
           Scalar m0,
           Scalar m1)
    { this->set(nSamples, points, m0, m1); }

    // convenience constructor
    template <class ScalarContainer>
    Spline(const ScalarContainer &x,
           const ScalarContainer &y,
           Scalar m0,
           Scalar m1)
    { this->set(x, y, m0, m1); }

    // convenience constructor
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

    // convenience constructor
    template <class ScalarContainer>
    Spline(const ScalarContainer &x,
           const ScalarContainer &y,
           Scalar m0, Scalar m1)
    {
        set(x,y,m0,m1);
    }
    
    // convenience constructor
    template <class XYContainer>
    Spline(const XYContainer &points,
           Scalar m0,
           Scalar m1)
    { this->set(points, m0, m1); }

    // convenience constructor
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
     * \brief Set the sampling points using an array of (x,y) tuples.
     *
     * Variant for a full spline.
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
