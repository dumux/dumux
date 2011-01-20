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
 * \brief Implements a spline with a variable number of sampling points
 */
#ifndef DUMUX_VARIABLE_LENGTH_SPLINE_HH
#define DUMUX_VARIABLE_LENGTH_SPLINE_HH

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/btdmatrix.hh>

#include "splinecommon_.hh"

namespace Dumux
{
//! \cond INTERNAL
/*
 * \brief The common code for all 3rd order polynomial splines with
 *        where the number of sampling points only known at run-time.
 */
template<class ScalarT>
class VariableLengthSpline_
    : public SplineCommon_<ScalarT,
                           VariableLengthSpline_<ScalarT> >
{
    friend class SplineCommon_<ScalarT, VariableLengthSpline_<ScalarT> >;

    typedef ScalarT Scalar;
    typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > Vector;
    typedef Dune::BTDMatrix<Dune::FieldMatrix<Scalar, 1, 1> > BTDMatrix;

public:
    /*!
     * \brief Returns the number of sampling points.
     */
    int numSamples() const
    { return xPos_.size(); }

    /*!
     * \brief Set the sampling points and the boundary slopes of the
     *        spline.
     */
    template <class ScalarContainer>
    void set(int nSamples,
             const ScalarContainer &x,
             const ScalarContainer &y,
             Scalar m0, Scalar m1)
    {
        xPos_.resize(nSamples);
        yPos_.resize(nSamples);
        m_.resize(nSamples);
        BTDMatrix M(nSamples);
        Vector d(nSamples);

        this->assignSamplingPoints_(xPos_, yPos_, x, y, numSamples());
        this->makeFullSystem_(M, d, m0, m1);

        // solve for the moments
        M.solve(m_, d);
    }

    /*!
     * \brief Set the sampling points and the boundary slopes of the
     *        spline.
     *
     * This assumes that the ScalarContainer type has a size method.
     */
    template <class ScalarContainer>
    void set(const ScalarContainer &x,
             const ScalarContainer &y,
             Scalar m0, Scalar m1)
    {
        assert(x.size() == y.size());
        set(x.size(), x, y, m0, m1);
    }

    /*!
     * \brief Set the sampling points for a full spline using a container of (x,y) tuples.
     */
    template <class XYContainer>
    void set(int nSamples,
             const XYContainer &xy,
             Scalar m0, Scalar m1)
    {
        xPos_.resize(nSamples);
        yPos_.resize(nSamples);
        m_.resize(nSamples);
        BTDMatrix M(nSamples);
        Vector d(nSamples);

        this->assignSamplingPoints_(xPos_, yPos_, xy, numSamples());
        this->makeFullSystem_(M, d, m0, m1);

        // solve for the moments
        M.solve(m_, d);
    }

    /*!
     * \brief Set the sampling points for a full spline using a container of (x,y) tuples.
     */
    template <class XYContainer>
    void set(const XYContainer &xy,
             Scalar m0, Scalar m1)
    {
        set(xy.size(), xy, m0, m1);
    }

    /*!
     * \brief Set the sampling points for a natural spline.
     */
    template <class ScalarContainer>
    void set(int nSamples,
             const ScalarContainer &x,
             const ScalarContainer &y)
    {
        xPos_.resize(nSamples);
        yPos_.resize(nSamples);
        m_.resize(nSamples);
        BTDMatrix M(nSamples);
        Vector d(nSamples);

        this->assignSamplingPoints_(xPos_, yPos_, x, y, numSamples());
        this->makeNaturalSystem_(M, d);

        // solve for the moments
        M.solve(m_, d);
    }

    /*!
     * \brief Set the sampling points for a natural spline.
     *
     * This assimes the ScalarContainer type has a size() method
     */
    template <class ScalarContainer>
    void set(const ScalarContainer &x,
             const ScalarContainer &y)
    {
        assert(x.size() == y.size());
        set(x.size(), x, y);
    }

    /*!
     * \brief Set the sampling points for a natural spline using a container of (x,y) tuples.
     */
    template <class XYContainer>
    void set(int nSamples,
             const XYContainer &points)
    {
        xPos_.resize(nSamples);
        yPos_.resize(nSamples);
        m_.resize(nSamples);
        BTDMatrix M(nSamples);
        Vector d(nSamples);

        this->assignSamplingPoints_(xPos_, yPos_, points, numSamples());
        this->makeNaturalSystem_(M, d);

        // solve for the moments
        M.solve(m_, d);
    }

    /*!
     * \brief Set the sampling points for a natural spline using a container of (x,y) tuples.
     *
     * This assimes the XYContainer type has a size() method
     */
    template <class XYContainer>
    void set(const XYContainer &points)
    { set(points.size(), points); }

protected:
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
//! \endcond
}

#endif
