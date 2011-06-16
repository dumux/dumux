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
 * \brief Implements a spline with a fixed number of sampling points
 */
#ifndef DUMUX_FIXED_LENGTH_SPLINE_HH
#define DUMUX_FIXED_LENGTH_SPLINE_HH

#include <dune/common/fvector.hh>
#include <dune/istl/btdmatrix.hh>

#include "splinecommon_.hh"

namespace Dumux
{
//! \cond INTERNAL
/*!
 * \brief The common code for all 3rd order polynomial splines with
 *        more than two sampling points.
 */
template<class ScalarT, int nSamples>
class FixedLengthSpline_
    : public SplineCommon_<ScalarT,
                           FixedLengthSpline_<ScalarT, nSamples> >
{
    friend class SplineCommon_<ScalarT, FixedLengthSpline_<ScalarT, nSamples> >;

    typedef ScalarT Scalar;
    typedef Dune::FieldVector<Scalar, nSamples> Vector;
    typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > BlockVector;
    typedef Dune::BTDMatrix<Dune::FieldMatrix<Scalar, 1, 1> > BTDMatrix;

protected:
    FixedLengthSpline_()
        : m_(nSamples)
    {};

public:
    /*!
     * \brief Returns the number of sampling points.
     */
    int numSamples() const
    { return nSamples; }

    /*!
     * \brief Set the sampling points and the boundary slopes of the
     *        spline.
     */
    template <class ScalarContainer>
    void set(const ScalarContainer &x,
             const ScalarContainer &y,
             Scalar m0, Scalar m1)
    {
        this->assignSamplingPoints_(xPos_, yPos_, x, y, numSamples());
        myMakeFullSpline_(m0, m1);
    }

    /*!
     * \brief Set the sampling points using an array of (x,y) pairs.
     *
     * Variant for a full spline.
     */
    template <class XYContainer>
    void setViaArrayContainer(const XYContainer &points, Scalar m0, Scalar m1)
    {
        assert(points.size() == numSamples());

        this->assignFromArrayList_(xPos_,
                                   yPos_,
                                   points.begin(),
                                   points.end(),
                                   numSamples());
        myMakeFullSpline_(m0, m1);
    }

    /*!
     * \brief Set the sampling points using an array of (x,y) pairs.
     *
     * Variant for a full spline.
     */
    template <class XYContainer>
    void setViaArrayArray(int nuSamples, const XYContainer &points, Scalar m0, Scalar m1)
    {
        assert(nuSamples == numSamples());

        this->assignFromArrayList_(xPos_,
                                   yPos_,
                                   &points[0],
                                   &points[numSamples() - 1],
                                   numSamples());
        myMakeFullSpline_(m0, m1);
    }

    /*!
     * \brief Set the sampling points using an array of (x,y) tuples.
     *
     * Variant for a full spline.
     */
    template <class TupleContainer>
    void setViaTupleContainer(const TupleContainer &points, Scalar m0, Scalar m1)
    {
        this->assignFromTupleList_(xPos_,
                                   yPos_,
                                   points.begin(),
                                   points.end(),
                                   numSamples());
        myMakeFullSpline_(m0, m1);
    }

    /*!
     * \brief Set the sampling points for a natural spline.
     */
    template <class ScalarContainer>
    void set(const ScalarContainer &x, const ScalarContainer &y)
    {
        this->assignSamplingPoints_(xPos_, yPos_, x, y, numSamples());
        myMakeNaturalSpline_();
    }

    /*!
     * \brief Set the sampling points using an array of (x,y) tuples.
     *
     * Variant for a natural spline.
     */
    template <class XYContainer>
    void setViaArrayContainer(const XYContainer &points)
    {
        this->assignFromArrayList_(xPos_,
                                   yPos_,
                                   points.begin(),
                                   points.end(),
                                   numSamples());
        myMakeNaturalSpline_();
    }

    /*!
     * \brief Set the sampling points using an array of (x,y) tuples.
     *
     * Variant for a natural spline.
     */
    template <class XYContainer>
    void setViaArrayArray(int nuSamples, const XYContainer &points)
    {
        assert(nuSamples == numSamples());

        this->assignFromArrayList_(xPos_,
                                   yPos_,
                                   &points[0],
                                   &points[numSamples() - 1], 
                                   numSamples());
        myMakeNaturalSpline_();
    }

    /*!
     * \brief Set the sampling points using an array of (x,y) tuples.
     *
     * Variant for a natural spline.
     */
    template <class TupleContainer>
    void setViaTupleContainer(const TupleContainer &points)
    {
        this->assignFromTupleList_(xPos_,
                                   yPos_,
                                   points.begin(),
                                   points.end(),
                                   numSamples());
        myMakeNaturalSpline_();
    }

protected:
    /*!
     * \brief Create a full spline from the already set sampling points.
     *
     * Also creates temporary matrix and right hand side vector.
     */
    void myMakeFullSpline_(Scalar m0, Scalar m1)
    {
        BTDMatrix M(numSamples());
        BlockVector d(numSamples());

        // create linear system of equations
        this->makeFullSystem_(M, d, m0, m1);

        // solve for the moments
        M.solve(m_, d);
    }

    /*!
     * \brief Create a natural spline from the already set sampling points.
     *
     * Also creates temporary matrix and right hand side vector.
     */
    void myMakeNaturalSpline_()
    {
        BTDMatrix M(numSamples());
        BlockVector d(numSamples());

        // create linear system of equations
        this->makeNaturalSystem_(M, d);

        // solve for the moments
        M.solve(m_, d);
    }

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
    BlockVector m_;
};

//! \endcond

}

#endif
