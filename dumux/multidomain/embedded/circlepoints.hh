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

namespace Dumux {
namespace EmbeddedCoupling {

/*!
 * \ingroup EmbeddedCoupling
 * \brief returns a vector of points on a circle
 * \param points a vector of points to be filled
 * \param center the circle center
 * \param normal the normal to the circle plane
 * \param radius the circle radius
 * \param numPoints the number of points
 */
template<class GlobalPosition, class Scalar>
void circlePoints(std::vector<GlobalPosition>& points,
                  const std::vector<Scalar>& sincos,
                  const GlobalPosition& center,
                  const GlobalPosition& normal,
                  const Scalar radius,
                  const std::size_t numPoints = 20)
{
    using std::abs;
    using ctype = typename GlobalPosition::value_type;

    constexpr ctype eps = 1.5e-7;
    static_assert(GlobalPosition::dimension == 3, "Only implemented for world dimension 3");

    // resize the points vector
    points.resize(numPoints);

    // make sure n is a unit vector
    auto n = normal;
    n /= n.two_norm();

    // caculate a vector u perpendicular to n
    GlobalPosition u;
    if (abs(n[0]) < eps && abs(n[1]) < eps)
        if (abs(n[2]) < eps)
            DUNE_THROW(Dune::MathError, "The normal vector has to be non-zero!");
        else
            u = {0, 1, 0};
    else
        u = {-n[1], n[0], 0};

    u *= radius/u.two_norm();

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
    using std::sin; using std::cos;
    std::vector<Scalar> sincos(2*numPoints);
    Scalar t = 0 + 0.1; // start with an offset
    for (std::size_t i = 0; i < numPoints; ++i)
    {
        sincos[2*i] = sin(t);
        sincos[2*i + 1] = cos(t);
        t += 2*M_PI/numPoints;
        if(t > 2*M_PI) t -= 2*M_PI;
    }

    circlePoints(points, sincos, center, normal, radius, numPoints);
    return points;
}

/*!
 * \ingroup EmbeddedCoupling
 * \brief returns a vector of points on a circle
 * \param points a vector of points to be filled
 * \param center the circle center
 * \param normal the normal to the circle plane
 * \param radius the circle radius
 * \param numPoints the number of points
 */
template<class GlobalPosition, class Scalar>
void ellipsePoints(std::vector<GlobalPosition>& points,
                   const std::vector<Scalar>& sincos,
                   const GlobalPosition& center,
                   const GlobalPosition& normal,
                   const Scalar a, const Scalar b,
                   const std::size_t numPoints = 20)
{
    using std::abs;
    using ctype = typename GlobalPosition::value_type;

    constexpr ctype eps = 1.5e-7;
    static_assert(GlobalPosition::dimension == 3, "Only implemented for world dimension 3");

    // resize the points vector
    points.resize(numPoints);

    // make sure n is a unit vector
    auto n = normal;
    n /= n.two_norm();

    // caculate a vector u perpendicular to n
    GlobalPosition u;
    if (abs(n[0]) < eps && abs(n[1]) < eps)
        if (abs(n[2]) < eps)
            DUNE_THROW(Dune::MathError, "The normal vector has to be non-zero!");
        else
            u = {0, 1, 0};
    else
        u = {-n[1], n[0], 0};

    u /= u.two_norm();

    // the ellipse parameterization is p(t) = b*cos(t)*u + a*sin(t)*(n x u) + c
    auto tangent = crossProduct(u, n);
    tangent /= tangent.two_norm();

    // insert the vertices
    for (std::size_t i = 0; i < numPoints; ++i)
    {
        points[i] = GlobalPosition({b*u[0]*sincos[2*i+1] + a*tangent[0]*sincos[2*i] + center[0],
                                    b*u[1]*sincos[2*i+1] + a*tangent[1]*sincos[2*i] + center[1],
                                    b*u[2]*sincos[2*i+1] + a*tangent[2]*sincos[2*i] + center[2]});
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
std::vector<GlobalPosition> ellipsePoints(const GlobalPosition& center,
                                          const GlobalPosition& normal,
                                          const Scalar a, const Scalar b,
                                          const std::size_t numPoints = 20)
{
    std::vector<GlobalPosition> points;

    // precompute the sin/cos
    using std::sin; using std::cos;
    std::vector<Scalar> sincos(2*numPoints);
    Scalar t = 0 + 0.1; // start with an offset
    for (std::size_t i = 0; i < numPoints; ++i)
    {
        sincos[2*i] = sin(t);
        sincos[2*i + 1] = cos(t);
        t += 2*M_PI/numPoints;
        if(t > 2*M_PI) t -= 2*M_PI;
    }

    ellipsePoints(points, sincos, center, normal, a, b, numPoints);
    return points;
}

} // end namespace EmbeddedCoupling
} // end namespace Dumux

#endif
