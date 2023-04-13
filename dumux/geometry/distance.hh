// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Geometry
 * \brief Helper functions for distance queries
 */
#ifndef DUMUX_GEOMETRY_DISTANCE_HH
#define DUMUX_GEOMETRY_DISTANCE_HH

#include <dune/common/fvector.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dumux/common/math.hh>
#include <dumux/geometry/boundingboxtree.hh>

namespace Dumux {

/*!
 * \ingroup Geometry
 * \brief Compute the average distance from a point to a geometry by integration
 */
template<class Geometry>
static inline typename Geometry::ctype
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
 * \brief Compute the squared distance from a point to a line through the points a and b
 */
template<class Point>
static inline typename Point::value_type
squaredDistancePointLine(const Point& p, const Point& a, const Point& b)
{
    const auto ab = b - a;
    const auto t = (p - a)*ab/ab.two_norm2();
    auto proj = a;
    proj.axpy(t, ab);
    return (proj - p).two_norm2();
}

/*!
 * \ingroup Geometry
 * \brief Compute the squared distance from a point to a line given by a geometry with two corners
 * \note We currently lack the a representation of a line geometry. This convenience function
 *       assumes a segment geometry (with two corners) is passed which represents a line geometry.
 */
template<class Geometry>
static inline typename Geometry::ctype
squaredDistancePointLine(const typename Geometry::GlobalCoordinate& p, const Geometry& geometry)
{
    static_assert(Geometry::mydimension == 1, "Geometry has to be a line");
    const auto& a = geometry.corner(0);
    const auto& b = geometry.corner(1);
    return squaredDistancePointLine(p, a, b);
}

/*!
 * \ingroup Geometry
 * \brief Compute the distance from a point to a line through the points a and b
 */
template<class Point>
static inline typename Point::value_type
distancePointLine(const Point& p, const Point& a, const Point& b)
{ using std::sqrt; return sqrt(squaredDistancePointLine(p, a, b)); }

/*!
 * \ingroup Geometry
 * \brief Compute the distance from a point to a line given by a geometry with two corners
 * \note We currently lack the a representation of a line geometry. This convenience function
 *       assumes a segment geometry (with two corners) is passed which represents a line geometry.
 */
template<class Geometry>
static inline typename Geometry::ctype
distancePointLine(const typename Geometry::GlobalCoordinate& p, const Geometry& geometry)
{ using std::sqrt; return sqrt(squaredDistancePointLine(p, geometry)); }

/*!
 * \ingroup Geometry
 * \brief Compute the squared distance from a point to the segment connecting the points a and b
 */
template<class Point>
static inline typename Point::value_type
squaredDistancePointSegment(const Point& p, const Point& a, const Point& b)
{
    const auto ab = b - a;
    const auto ap = p - a;
    const auto t = ap*ab;

    if (t <= 0.0)
        return ap.two_norm2();

    const auto lengthSq = ab.two_norm2();
    if (t >= lengthSq)
        return (b - p).two_norm2();

    auto proj = a;
    proj.axpy(t/lengthSq, ab);
    return (proj - p).two_norm2();
}

/*!
 * \ingroup Geometry
 * \brief Compute the squared distance from a point to a given segment geometry
 */
template<class Geometry>
static inline typename Geometry::ctype
squaredDistancePointSegment(const typename Geometry::GlobalCoordinate& p, const Geometry& geometry)
{
    static_assert(Geometry::mydimension == 1, "Geometry has to be a segment");
    const auto& a = geometry.corner(0);
    const auto& b = geometry.corner(1);
    return squaredDistancePointSegment(p, a, b);
}

/*!
 * \ingroup Geometry
 * \brief Compute the distance from a point to the segment connecting the points a and b
 */
template<class Point>
static inline typename Point::value_type
distancePointSegment(const Point& p, const Point& a, const Point& b)
{ using std::sqrt; return sqrt(squaredDistancePointSegment(p, a, b)); }

/*!
 * \ingroup Geometry
 * \brief Compute the distance from a point to a given segment geometry
 */
template<class Geometry>
static inline typename Geometry::ctype
distancePointSegment(const typename Geometry::GlobalCoordinate& p, const Geometry& geometry)
{ using std::sqrt; return sqrt(squaredDistancePointSegment(p, geometry)); }

/*!
 * \ingroup Geometry
 * \brief Compute the shortest squared distance from a point to the triangle connecting the points a, b and c
 *        See https://www.iquilezles.org/www/articles/triangledistance/triangledistance.htm.
 */
template<class Point>
static inline typename Point::value_type
squaredDistancePointTriangle(const Point& p, const Point& a, const Point& b, const Point& c)
{
    static_assert(Point::dimension == 3, "Only works in 3D");
    const auto ab = b - a;
    const auto bc = c - b;
    const auto ca = a - c;
    const auto normal = crossProduct(ab, ca);

    const auto ap = p - a;
    const auto bp = p - b;
    const auto cp = p - c;

    const auto sum = sign(crossProduct(ab, normal)*ap)
                   + sign(crossProduct(bc, normal)*bp)
                   + sign(crossProduct(ca, normal)*cp);

    // if there is no orthogonal projection
    // (point is outside the infinite prism implied by the triangle and its surface normal)
    // compute distance to the edges (codim-1 facets)
    if (sum < 2.0)
    {
        using std::min;
        return min({squaredDistancePointSegment(p, a, b),
                    squaredDistancePointSegment(p, a, c),
                    squaredDistancePointSegment(p, b, c)});
    }
    // compute distance via orthogonal projection
    else
    {
        const auto tmp = normal*ap;
        return tmp*tmp / (normal*normal);
    }
}

/*!
 * \ingroup Geometry
 * \brief Compute the shortest squared distance from a point to a given triangle geometry
 */
template<class Geometry>
static inline typename Geometry::ctype
squaredDistancePointTriangle(const typename Geometry::GlobalCoordinate& p, const Geometry& geometry)
{
    static_assert(Geometry::coorddimension == 3, "Only works in 3D");
    static_assert(Geometry::mydimension == 2, "Geometry has to be a triangle");
    assert(geometry.corners() == 3);
    const auto& a = geometry.corner(0);
    const auto& b = geometry.corner(1);
    const auto& c = geometry.corner(2);
    return squaredDistancePointTriangle(p, a, b, c);
}

/*!
 * \ingroup Geometry
 * \brief Compute the shortest distance from a point to the triangle connecting the points a, b and c
 */
template<class Point>
static inline typename Point::value_type
distancePointTriangle(const Point& p, const Point& a, const Point& b, const Point& c)
{ using std::sqrt; return sqrt(squaredDistancePointTriangle(p, a, b, c)); }

/*!
 * \ingroup Geometry
 * \brief Compute the shortest distance from a point to a given triangle geometry
 */
template<class Geometry>
static inline typename Geometry::ctype
distancePointTriangle(const typename Geometry::GlobalCoordinate& p, const Geometry& geometry)
{ using std::sqrt; return sqrt(squaredDistancePointTriangle(p, geometry)); }

/*!
 * \ingroup Geometry
 * \brief Compute the shortest squared distance from a point to a given polygon geometry
 * \note We only support triangles and quadrilaterals so far.
 */
template<class Geometry>
static inline typename Geometry::ctype
squaredDistancePointPolygon(const typename Geometry::GlobalCoordinate& p, const Geometry& geometry)
{
    if (geometry.corners() == 3)
        return squaredDistancePointTriangle(p, geometry);
    else if (geometry.corners() == 4)
    {
        const auto& a = geometry.corner(0);
        const auto& b = geometry.corner(1);
        const auto& c = geometry.corner(2);
        const auto& d = geometry.corner(3);

        using std::min;
        return min(squaredDistancePointTriangle(p, a, b, d),
                   squaredDistancePointTriangle(p, a, d, c));
    }
    else
        DUNE_THROW(Dune::NotImplemented, "Polygon with " << geometry.corners() << " corners not supported");
}

/*!
 * \ingroup Geometry
 * \brief Compute the shortest distance from a point to a given polygon geometry
 * \note We only support triangles and quadrilaterals so far.
 */
template<class Geometry>
static inline typename Geometry::ctype
distancePointPolygon(const typename Geometry::GlobalCoordinate& p, const Geometry& geometry)
{ using std::sqrt; return sqrt(squaredDistancePointPolygon(p, geometry)); }

/*!
 * \ingroup Geometry
 * \brief Compute the average distance from a segment to a geometry by integration
 */
template<class Geometry>
static inline typename Geometry::ctype
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
static inline ctype distance(const Dune::FieldVector<ctype, dimWorld>& a,
                             const Dune::FieldVector<ctype, dimWorld>& b)
{ return (a-b).two_norm(); }

/*!
 * \ingroup Geometry
 * \brief Compute the shortest squared distance between two points
 */
template<class ctype, int dimWorld>
static inline ctype squaredDistance(const Dune::FieldVector<ctype, dimWorld>& a,
                                    const Dune::FieldVector<ctype, dimWorld>& b)
{ return (a-b).two_norm2(); }


namespace Detail {

// helper struct to compute distance between two geometries, specialized below
template<class Geo1, class Geo2,
         int dimWorld = Geo1::coorddimension,
         int dim1 = Geo1::mydimension, int dim2 = Geo2::mydimension>
struct GeometrySquaredDistance
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
struct GeometrySquaredDistance<Geo1, Geo2, dimWorld, 0, 0>
{
    static_assert(Geo1::coorddimension == Geo2::coorddimension, "Geometries have to have the same coordinate dimensions");
    static auto distance(const Geo1& geo1, const Geo2& geo2)
    { return Dumux::squaredDistance(geo1.corner(0), geo2.corner(0)); }
};

// distance segment-point
template<class Geo1, class Geo2, int dimWorld>
struct GeometrySquaredDistance<Geo1, Geo2, dimWorld, 1, 0>
{
    static_assert(Geo1::coorddimension == Geo2::coorddimension, "Geometries have to have the same coordinate dimensions");
    static auto distance(const Geo1& geo1, const Geo2& geo2)
    { return squaredDistancePointSegment(geo2.corner(0), geo1); }
};

// distance point-segment
template<class Geo1, class Geo2, int dimWorld>
struct GeometrySquaredDistance<Geo1, Geo2, dimWorld, 0, 1>
{
    static_assert(Geo1::coorddimension == Geo2::coorddimension, "Geometries have to have the same coordinate dimensions");
    static auto distance(const Geo1& geo1, const Geo2& geo2)
    { return squaredDistancePointSegment(geo1.corner(0), geo2); }
};

// distance point-polygon
template<class Geo1, class Geo2, int dimWorld>
struct GeometrySquaredDistance<Geo1, Geo2, dimWorld, 0, 2>
{
    static_assert(Geo1::coorddimension == Geo2::coorddimension, "Geometries have to have the same coordinate dimensions");
    static inline auto distance(const Geo1& geo1, const Geo2& geo2)
    { return squaredDistancePointPolygon(geo1.corner(0), geo2); }
};

// distance polygon-point
template<class Geo1, class Geo2, int dimWorld>
struct GeometrySquaredDistance<Geo1, Geo2, dimWorld, 2, 0>
{
    static_assert(Geo1::coorddimension == Geo2::coorddimension, "Geometries have to have the same coordinate dimensions");
    static inline auto distance(const Geo1& geo1, const Geo2& geo2)
    { return squaredDistancePointPolygon(geo2.corner(0), geo1); }
};

} // end namespace Detail

/*!
 * \ingroup Geometry
 * \brief Compute the shortest distance between two geometries
 */
template<class Geo1, class Geo2>
static inline auto squaredDistance(const Geo1& geo1, const Geo2& geo2)
{ return Detail::GeometrySquaredDistance<Geo1, Geo2>::distance(geo1, geo2); }

/*!
 * \ingroup Geometry
 * \brief Compute the shortest distance between two geometries
 */
template<class Geo1, class Geo2>
static inline auto distance(const Geo1& geo1, const Geo2& geo2)
{ using std::sqrt; return sqrt(squaredDistance(geo1, geo2)); }


//! Distance implementation details
namespace Detail {

/*!
 * \ingroup Geometry
 * \brief Compute the closest entity in an AABB tree (index and shortest squared distance) recursively
 * \note specialization for geometries with dimension larger than points
 */
template<class EntitySet, class ctype, int dimworld,
         typename std::enable_if_t<(EntitySet::Entity::Geometry::mydimension > 0), int> = 0>
void closestEntity(const Dune::FieldVector<ctype, dimworld>& point,
                   const BoundingBoxTree<EntitySet>& tree,
                   std::size_t node,
                   ctype& minSquaredDistance,
                   std::size_t& eIdx)
{
    // Get the bounding box for the current node
    const auto& bBox = tree.getBoundingBoxNode(node);

    // If bounding box is outside radius, then don't search further
    const auto squaredDistance = squaredDistancePointBoundingBox(
        point, tree.getBoundingBoxCoordinates(node)
    );

    // do not continue if the AABB is further away than the current minimum
    if (squaredDistance > minSquaredDistance) return;

    // If the bounding box is a leaf, update distance and index with primitive test
    else if (tree.isLeaf(bBox, node))
    {
        const std::size_t entityIdx = bBox.child1;
        const auto squaredDistance = [&]{
            const auto geometry = tree.entitySet().entity(entityIdx).geometry();
            if constexpr (EntitySet::Entity::Geometry::mydimension == 2)
                return squaredDistancePointPolygon(point, geometry);
            else if constexpr (EntitySet::Entity::Geometry::mydimension == 1)
                return squaredDistancePointSegment(point, geometry);
            else
                DUNE_THROW(Dune::NotImplemented, "squaredDistance to entity with dim>2");
        }();

        if (squaredDistance < minSquaredDistance)
        {
            eIdx = entityIdx;
            minSquaredDistance = squaredDistance;
        }
    }

    // Check the children nodes recursively
    else
    {
        closestEntity(point, tree, bBox.child0, minSquaredDistance, eIdx);
        closestEntity(point, tree, bBox.child1, minSquaredDistance, eIdx);
    }
}

/*!
 * \ingroup Geometry
 * \brief Compute the closest entity in an AABB tree (index and shortest squared distance) recursively
 * \note specialization for point geometries (point cloud AABB tree)
 */
template<class EntitySet, class ctype, int dimworld,
         typename std::enable_if_t<(EntitySet::Entity::Geometry::mydimension == 0), int> = 0>
void closestEntity(const Dune::FieldVector<ctype, dimworld>& point,
                   const BoundingBoxTree<EntitySet>& tree,
                   std::size_t node,
                   ctype& minSquaredDistance,
                   std::size_t& eIdx)
{
    // Get the bounding box for the current node
    const auto& bBox = tree.getBoundingBoxNode(node);

    // If the bounding box is a leaf, update distance and index with primitive test
    if (tree.isLeaf(bBox, node))
    {
        const std::size_t entityIdx = bBox.child1;
        const auto& p = tree.entitySet().entity(entityIdx).geometry().corner(0);
        const auto squaredDistance = (p-point).two_norm2();

        if (squaredDistance < minSquaredDistance)
        {
            eIdx = entityIdx;
            minSquaredDistance = squaredDistance;
        }
    }

    // Check the children nodes recursively
    else
    {
        // If bounding box is outside radius, then don't search further
        const auto squaredDistance = squaredDistancePointBoundingBox(
            point, tree.getBoundingBoxCoordinates(node)
        );

        // do not continue if the AABB is further away than the current minimum
        if (squaredDistance > minSquaredDistance) return;

        closestEntity(point, tree, bBox.child0, minSquaredDistance, eIdx);
        closestEntity(point, tree, bBox.child1, minSquaredDistance, eIdx);
    }
}

} // end namespace Detail

/*!
 * \ingroup Geometry
 * \brief Compute the closest entity in an AABB tree to a point (index and shortest squared distance)
 * \param point the point
 * \param tree the AABB tree
 * \param minSquaredDistance conservative estimate of the minimum distance
 *
 * \note it is important that if an estimate is provided for minSquaredDistance to choose
 *       this estimate to be larger than the expected result. If the minSquaredDistance is smaller
 *       or equal to the result the returned entity index will be zero and the distance is equal
 *       to the estimate. However, this can also be the correct result. When in doubt, use the
 *       default parameter value.
 */
template<class EntitySet, class ctype, int dimworld>
std::pair<ctype, std::size_t> closestEntity(const Dune::FieldVector<ctype, dimworld>& point,
                                            const BoundingBoxTree<EntitySet>& tree,
                                            ctype minSquaredDistance = std::numeric_limits<ctype>::max())
{
    std::size_t eIdx = 0;
    Detail::closestEntity(point, tree, tree.numBoundingBoxes() - 1, minSquaredDistance, eIdx);
    using std::sqrt;
    return { minSquaredDistance, eIdx };
}

/*!
 * \ingroup Geometry
 * \brief Compute the shortest squared distance to entities in an AABB tree
 */
template<class EntitySet, class ctype, int dimworld>
ctype squaredDistance(const Dune::FieldVector<ctype, dimworld>& point,
                      const BoundingBoxTree<EntitySet>& tree,
                      ctype minSquaredDistance = std::numeric_limits<ctype>::max())
{
    return closestEntity(point, tree, minSquaredDistance).first;
}

/*!
 * \ingroup Geometry
 * \brief Compute the shortest distance to entities in an AABB tree
 */
template<class EntitySet, class ctype, int dimworld>
ctype distance(const Dune::FieldVector<ctype, dimworld>& point,
               const BoundingBoxTree<EntitySet>& tree,
               ctype minSquaredDistance = std::numeric_limits<ctype>::max())
{
    using std::sqrt;
    return sqrt(squaredDistance(point, tree, minSquaredDistance));
}

} // end namespace Dumux

#endif
