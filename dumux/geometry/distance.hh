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
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Geometry
 * \brief Helper functions for distance queries
 */
#ifndef DUMUX_GEOMETRY_DISTANCE_HH
#define DUMUX_GEOMETRY_DISTANCE_HH

#include <dune/common/fvector.hh>
#include <dune/geometry/quadraturerules.hh>

namespace Dumux {

/*!
 * \ingroup Geometry
 * \brief Compute the average distance from a point to a geometry by integration
 */
template<class Geometry>
inline typename Geometry::ctype
averageDistancePointGeometry(const typename Geometry::GlobalCoordinate& p,
                             const Geometry& geometry,
                             std::size_t integrationOrder = 2)
{
    typename Geometry::ctype avgDist = 0.0;
    const auto& quad = Dune::QuadratureRules<typename Geometry::ctype, Geometry::mydimension>::rule(geometry.type(), integrationOrder);
    for (const auto& qp : quad)
        avgDist += (geometry.global(qp.position())-p).two_norm()*qp.weight()*geometry.integrationElement(qp.position());
    return avgDist/geometry.volume();
}

/*!
 * \ingroup Geometry
 * \brief Compute the distance from a point to a line through the points a and b
 */
template<class Point>
inline typename Point::value_type
distancePointLine(const Point& p, const Point& a, const Point& b)
{
    const auto ab = b - a;
    const auto t = (p - a)*ab/ab.two_norm2();
    auto proj = a;
    proj.axpy(t, ab);
    return (proj - p).two_norm();
}

/*!
 * \ingroup Geometry
 * \brief Compute the distance from a point to a line given by a geometry with two corners
 * \note We currently lack the a representation of a line geometry. This convenience function
 *       assumes a segment geometry (with two corners) is passed which represents a line geometry.
 */
template<class Geometry>
inline typename Geometry::ctype
distancePointLine(const typename Geometry::GlobalCoordinate& p, const Geometry& geometry)
{
    static_assert(Geometry::mydimension == 1, "Geometry has to be a line");
    const auto& a = geometry.corner(0);
    const auto& b = geometry.corner(1);
    return distancePointLine(p, a, b);
}

/*!
 * \ingroup Geometry
 * \brief Compute the distance from a point to the segment connecting the points a and b
 */
template<class Point>
inline typename Point::value_type
distancePointSegment(const Point& p, const Point& a, const Point& b)
{
    const auto ab = b - a;
    auto t = (p - a)*ab;

    if (t <= 0.0)
        return (a - p).two_norm();

    const auto lengthSq = ab.two_norm2();
    if (t >= lengthSq)
        return (b - p).two_norm();

    auto proj = a;
    proj.axpy(t/lengthSq, ab);
    return (proj - p).two_norm();
}

/*!
 * \ingroup Geometry
 * \brief Compute the distance from a point to a given segment geometry
 */
template<class Geometry>
inline typename Geometry::ctype
distancePointSegment(const typename Geometry::GlobalCoordinate& p, const Geometry& geometry)
{
    static_assert(Geometry::mydimension == 1, "Geometry has to be a segment");
    const auto& a = geometry.corner(0);
    const auto& b = geometry.corner(1);
    return distancePointSegment(p, a, b);
}

/*!
 * \ingroup Geometry
 * \brief Compute the average distance from a segment to a geometry by integration
 */
template<class Geometry>
inline typename Geometry::ctype
averageDistanceSegmentGeometry(const typename Geometry::GlobalCoordinate& a,
                               const typename Geometry::GlobalCoordinate& b,
                               const Geometry& geometry,
                               std::size_t integrationOrder = 2)
{
    typename Geometry::ctype avgDist = 0.0;
    const auto& quad = Dune::QuadratureRules<typename Geometry::ctype, Geometry::mydimension>::rule(geometry.type(), integrationOrder);
    for (const auto& qp : quad)
        avgDist += distancePointSegment(geometry.global(qp.position()), a, b)*qp.weight()*geometry.integrationElement(qp.position());
    return avgDist/geometry.volume();
}

/*!
 * \ingroup Geometry
 * \brief Compute the shortest distance between two points
 */
template<class ctype, int dimWorld>
inline ctype distance(const Dune::FieldVector<ctype, dimWorld>& a,
                      const Dune::FieldVector<ctype, dimWorld>& b)
{ return (a-b).two_norm(); }



namespace Detail {

// helper struct to compute distance between two geometries, specialized below
template<class Geo1, class Geo2,
         int dimWorld = Geo1::coorddimension,
         int dim1 = Geo1::mydimension, int dim2 = Geo2::mydimension>
struct GeometryDistance
{
    static_assert(Geo1::coorddimension == Geo2::coorddimension, "Geometries have to have the same coordinate dimensions");
    static auto distance(const Geo1& geo1, const Geo2& geo2)
    {
        DUNE_THROW(Dune::NotImplemented, "Geometry distance computation not implemented for dimworld = "
                    << dimWorld << ", dim1 = " << dim1 << ", dim2 = " << dim2);
    }
};

// distance point-point
template<class Geo1, class Geo2, int dimWorld>
struct GeometryDistance<Geo1, Geo2, dimWorld, 0, 0>
{
    static_assert(Geo1::coorddimension == Geo2::coorddimension, "Geometries have to have the same coordinate dimensions");
    static auto distance(const Geo1& geo1, const Geo2& geo2)
    { return Dumux::distance(geo1.corner(0), geo2.corner(0)); }
};

// distance segment-point
template<class Geo1, class Geo2, int dimWorld>
struct GeometryDistance<Geo1, Geo2, dimWorld, 1, 0>
{
    static_assert(Geo1::coorddimension == Geo2::coorddimension, "Geometries have to have the same coordinate dimensions");
    static auto distance(const Geo1& geo1, const Geo2& geo2)
    { return distancePointSegment(geo2.corner(0), geo1); }
};

// distance point-segment
template<class Geo1, class Geo2, int dimWorld>
struct GeometryDistance<Geo1, Geo2, dimWorld, 0, 1>
{
    static_assert(Geo1::coorddimension == Geo2::coorddimension, "Geometries have to have the same coordinate dimensions");
    static auto distance(const Geo1& geo1, const Geo2& geo2)
    { return distancePointSegment(geo1.corner(0), geo2); }
};

} // end namespace Detail

/*!
 * \ingroup Geometry
 * \brief Compute the shortest distance between two geometries
 */
template<class Geo1, class Geo2>
inline auto distance(const Geo1& geo1, const Geo2& geo2)
{ return Detail::GeometryDistance<Geo1, Geo2>::distance(geo1, geo2); }

} // end namespace Dumux

#endif
