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
 * \ingroup Geomechanics
 * \brief MohrCircle class
 */
#ifndef DUMUX_STRESS_STATE_MOHR_CIRCLE_SPACE_HH
#define DUMUX_STRESS_STATE_MOHR_CIRCLE_SPACE_HH

#include <cmath>
#include <array>
#include "point.hh"

namespace Dumux{
template<class Scalar, class Point, class StressTensor>
class MohrCircle{
using PointList = Dune::FieldVector<Point,2>;
public:
    /*!
     * \brief create Mohr Circle from stress tensor.
     *  sign convention see Figure 5 Option 3 in https://en.wikipedia.org/wiki/Mohr%27s_circle
     * \param stress
     */
    MohrCircle(const StressTensor& stress)
    : circleCenter_{(stress[0][0]+stress[1][1])/2.0,0.0}
    {
        Point p1(stress[0][0],-stress[0][1]);
        auto dist = p1.pos() - circleCenter_.pos();
        radius_ = dist.two_norm();
    }

    Scalar radius() const {return radius_;}
    Point center() const {return circleCenter_;}

    /*!
     * \brief return the intersections of line and the Mohr-circle
     *  Point[0].x >= Point[1].x
     * \tparam Line
     * \param l line
     * \return PointList
     */
    template<class Line>
    PointList intersections(const Line& l) const
    {
        // solve equation (x-x_c)^2 + (ax+b)^2 = r^2
        using std::sqrt;
        Scalar co_a = 1 + l.slope()*l.slope();
        Scalar co_b = 2 * (l.slope()* l.intercept() - circleCenter_.x());
        Scalar co_c = l.intercept() * l.intercept() + circleCenter_.x()*circleCenter_.x() - radius_*radius_;
        Scalar delta = sqrt(co_b*co_b - 4*co_a*co_c);
        Scalar x0 = (-co_b + delta )/(2*co_a);
        Scalar x1 = (-co_b - delta )/(2*co_a);
        Scalar y0 = l.y(x0);
        Scalar y1 = l.y(x1);

        PointList intersections;
        intersections[0] = Point{x0, y0};
        intersections[1] = Point{x1, y1};

        return intersections;
    }

    template<class Line>
    bool hasIntersection(const Line& l) const
    {
        Scalar dist = pointLineDistance<Scalar>(circleCenter_, l);
        return dist < radius_;
    }

private:
    Point circleCenter_;
    Scalar radius_;
};
}// end namespace Dumux
#endif
