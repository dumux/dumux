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
 * \ingroup EmbeddedCoupling
 * \brief Helper function to compute points on a circle
 */

#ifndef DUMUX_MULTIDOMAIN_EMBEDDED_CIRCLEPOINTS_HH
#define DUMUX_MULTIDOMAIN_EMBEDDED_CIRCLEPOINTS_HH

#include <vector>
#include <cmath>

#include <dune/common/exceptions.hh>

#include <dumux/common/math.hh>
#include <dumux/geometry/normal.hh>

namespace Dumux::EmbeddedCoupling {

/*!
 * \ingroup EmbeddedCoupling
 * \param points a vector of points to be filled
 * \param sincos vector with [sin(a0), cos(a0), sin(a1), cos(a1), ...] for each circumferential sample point ai
 * \param center the circle center
 * \param normal the normal to the circle plane
 * \param radius the circle radius
 * \note This version allows for a more efficient circle point generator
 *       if the circumferential positions are fixed and can be reused.
 *       In this case, sine and cosine of the corresponding angles can be precomputed.
 */
template<class GlobalPosition, class Scalar>
void circlePoints(std::vector<GlobalPosition>& points,
                  const std::vector<Scalar>& sincos,
                  const GlobalPosition& center,
                  const GlobalPosition& normal,
                  const Scalar radius)
{
    assert(sincos.size() % 2 == 0 && "Sample angles have to be pairs of sin/cos so size needs to be even.");
    static_assert(GlobalPosition::dimension == 3, "Only implemented for world dimension 3");

    // resize the points vector
    const std::size_t numPoints = sincos.size()/2;
    points.resize(numPoints);

    // make sure n is a unit vector
    auto n = normal;
    n /= n.two_norm();

    // calculate a vector u perpendicular to n
    auto u = unitNormal(n); u*= radius;

    // the circle parameterization is p(t) = r*cos(t)*u + r*sin(t)*(n x u) + c
    auto tangent = crossProduct(u, n);
    tangent *= radius/tangent.two_norm();

    // insert the vertices
    for (std::size_t i = 0; i < numPoints; ++i)
    {
        points[i] = GlobalPosition({u[0]*sincos[2*i+1] + tangent[0]*sincos[2*i] + center[0],
                                    u[1]*sincos[2*i+1] + tangent[1]*sincos[2*i] + center[1],
                                    u[2]*sincos[2*i+1] + tangent[2]*sincos[2*i] + center[2]});
    }
}

/*!
 * \ingroup EmbeddedCoupling
 * \brief returns a vector of points on a circle
 * \param center the circle center
 * \param normal the normal to the circle plane
 * \param radius the circle radius
 * \param numPoints the number of points
 */
template<class GlobalPosition, class Scalar>
std::vector<GlobalPosition> circlePoints(const GlobalPosition& center,
                                         const GlobalPosition& normal,
                                         const Scalar radius,
                                         const std::size_t numPoints = 20)
{
    std::vector<GlobalPosition> points;

    // precompute the sin/cos
    std::vector<Scalar> sincos(2*numPoints);
    Scalar t = 0 + 0.1; // add arbitrary offset
    for (std::size_t i = 0; i < numPoints; ++i)
    {
        using std::sin; using std::cos;
        sincos[2*i] = sin(t);
        sincos[2*i + 1] = cos(t);
        t += 2*M_PI/numPoints;
        if(t > 2*M_PI) t -= 2*M_PI;
    }

    circlePoints(points, sincos, center, normal, radius);
    return points;
}

} // end namespace Dumux::EmbeddedCoupling

#endif
