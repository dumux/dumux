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
 * \brief A class for collision detection of two geometries
 *        and computation of intersection corners
 */
#ifndef DUMUX_GEOMETRY_INTERSECTION_HH
#define DUMUX_GEOMETRY_INTERSECTION_HH

#include <dune/common/exceptions.hh>
#include <dune/common/promotiontraits.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/math.hh>
#include <dumux/common/geometry/intersectspointgeometry.hh>
#include <dumux/common/geometry/grahamconvexhull.hh>

namespace Dumux {

/*!
 * \ingroup Geometry
 * \brief A class for geometry collision detection and intersection calculation
 * The class can be specialized for combinations of dimworld, dim1, dim2, where
 * dimworld is the world dimension embedding a grid of dim1 and a grid of dim2.
 */
template
<class Geometry1, class Geometry2,
  int dimworld = Geometry1::coorddimension,
  int dim1 = Geometry1::mydimension,
  int dim2 = Geometry2::mydimension>
class GeometryIntersection
{
public:
    //! Determine if the two geometries intersect and compute the intersection corners
    template<class IntersectionType>
    static bool intersection(const Geometry1& geo1, const Geometry2& geo2, IntersectionType& intersection)
    {
        static_assert(dimworld == Geometry2::coorddimension, "Can only intersect geometries of same coordinate dimension");
        DUNE_THROW(Dune::NotImplemented, "Geometry intersection detection with intersection computation not implemented for dimworld = "
                                          << dimworld << ", dim1 = " << dim1 << ", dim2 = " << dim2);
    }
};

/*!
 * \ingroup Geometry
 * \brief A class for polygon--segment intersection in 2d space
 */
template <class Geometry1, class Geometry2>
class GeometryIntersection<Geometry1, Geometry2, 2, 2, 1>
{
    enum { dimworld = 2 };
    enum { dim1 = 2 };
    enum { dim2 = 1 };

public:
    using ctype = typename Dune::PromotionTraits<typename Geometry1::ctype, typename Geometry2::ctype>::PromotedType;
    using Point = Dune::FieldVector<ctype, dimworld>;
    using IntersectionType = std::array<std::vector<Point>, 1>;

private:
    static constexpr ctype eps_ = 1.5e-7; // base epsilon for floating point comparisons
    using ReferenceElements = typename Dune::ReferenceElements<ctype, dim1>;

public:
    /*!
     * \brief Colliding segment and convex polygon
     * \param geo1/geo2 The geometries to intersect
     * \param intersection If the geometries collide intersection holds the
     *        corner points of the intersection object in global coordinates.
     */
    static bool intersection(const Geometry1& geo1, const Geometry2& geo2, IntersectionType& intersection)
    {
        const auto a = geo2.corner(0);
        const auto b = geo2.corner(1);
        const auto d = b - a;

        const std::vector<std::array<int, 2>> edges = [&]()
        {
            std::vector<std::array<int, 2>> edges;
            switch (geo1.corners())
            {
                case 4: // quadrilateral
                    edges = {{0, 2}, {3, 1}, {1, 0}, {2, 3}};
                    break;
                case 3: // triangle
                    edges = {{1, 0}, {0, 2}, {2, 1}};
                    break;
                default:
                    DUNE_THROW(Dune::NotImplemented, "Collision of segment and geometry of type "
                                   << geo1.type() << " with "<< geo1.corners() << " corners.");
            }

            return edges;
        } ();

        // The initial interval is the whole segment
        // afterward we start clipping the interval
        // by the lines decribed by the edges
        ctype tfirst = 0.0;
        ctype tlast = 1.0;

        const auto center1 = geo1.center();
        for (const auto& e : edges)
        {
            // compute outer normal vector of the edge
            const auto c0 = geo1.corner(e[0]);
            const auto c1 = geo1.corner(e[1]);
            const auto edge = c1 - c0;
            const auto eps = eps_*c0.two_norm();

            Dune::FieldVector<ctype, dimworld> n({-1.0*edge[1], edge[0]});
            n /= n.two_norm();

            // orientation of the element is not known
            // make sure normal points outwards of element
            if ( n*(center1-c0) > 0.0 )
                n *= -1.0;

            const ctype denom = n*d;
            const ctype dist = n*(a-c0);

            // if denominator is zero the segment in parallel to the edge.
            // If the distance is positive there is no intersection
            using std::abs;
            if (abs(denom) < eps)
            {
                if (dist > eps)
                    return false;
            }
            else // not parallel: compute line-line intersection
            {
                const ctype t = -dist / denom;
                // if entering half space cut tfirst if t is larger
                using std::signbit;
                if (signbit(denom))
                {
                    if (t > tfirst)
                        tfirst = t;
                }
                // if exiting half space cut tlast if t is smaller
                else
                {
                    if (t < tlast)
                        tlast = t;
                }
                // there is no intersection of the interval is empty
                // use unscaled epsilon since t is in local coordinates
                if (tfirst > tlast - eps_)
                    return false;
            }
        }
        // If we made it until here an intersection exists. We also export
        // the intersections geometry now s(t) = a + t(b-a) in [tfirst, tlast]
        intersection = {std::vector<Point>({geo2.global(tfirst), geo2.global(tlast)})};
        return true;
    }
};

/*!
 * \ingroup Geometry
 * \brief A class for polyhedron--segment intersection in 3d space
 */
template <class Geometry1, class Geometry2>
class GeometryIntersection<Geometry1, Geometry2, 3, 3, 1>
{
    enum { dimworld = 3 };
    enum { dim1 = 3 };
    enum { dim2 = 1 };

public:
    using ctype = typename Dune::PromotionTraits<typename Geometry1::ctype, typename Geometry2::ctype>::PromotedType;
    using Point = Dune::FieldVector<ctype, dimworld>;
    using IntersectionType = std::array<std::vector<Point>, 1>;

private:
    static constexpr ctype eps_ = 1.5e-7; // base epsilon for floating point comparisons
    using ReferenceElements = typename Dune::ReferenceElements<ctype, dim1>;

public:
    /*!
     *  \brief Colliding segment and convex polyhedron
     *  \note Algorithm based on the one from "Real-Time Collision Detection" by Christer Ericson,
     *        published by Morgan Kaufmann Publishers, (c) 2005 Elsevier Inc.
     *        Basis is the theorem that for any two non-intersecting convex polyhedrons
     *        a separating plane exists.
     * \param geo1/geo2 The geometries to intersect
     * \param intersection If the geometries collide intersection holds the corner points of
     *        the intersection object in global coordinates.
     */
    static bool intersection(const Geometry1& geo1, const Geometry2& geo2, IntersectionType& intersection)
    {
        static_assert(int(dimworld) == int(Geometry2::coorddimension), "Can only collide geometries of same coordinate dimension");

        const auto a = geo2.corner(0);
        const auto b = geo2.corner(1);
        const auto d = b - a;

        // The initial interval is the whole segment
        // afterward we start clipping the interval
        // by the planes decribed by the facet
        ctype tfirst = 0.0;
        ctype tlast = 1.0;

        const std::vector<std::vector<int>> facets = [&]()
        {
            std::vector<std::vector<int>> facets;
            // sort facet corner so that normal n = (p1-p0)x(p2-p0) always points outwards
            switch (geo1.corners())
            {
                // todo compute intersection geometries
                case 8: // hexahedron
                    facets = {{2, 0, 6, 4}, {1, 3, 5, 7}, {0, 1, 4, 5},
                              {3, 2, 7, 6}, {1, 0, 3, 2}, {4, 5, 6, 7}};
                    break;
                case 4: // tetrahedron
                    facets = {{1, 0, 2}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3}};
                    break;
                default:
                    DUNE_THROW(Dune::NotImplemented, "Collision of segment and geometry of type "
                                   << geo1.type() << ", "<< geo1.corners() << " corners.");
            }

            return facets;
        }();

        for (const auto& f : facets)
        {
            // compute normal vector by cross product
            const auto v0 = geo1.corner(f[1]) - geo1.corner(f[0]);
            const auto v1 = geo1.corner(f[2]) - geo1.corner(f[0]);
            const auto eps = eps_*v0.two_norm();

            auto n = crossProduct(v0, v1);
            n /= n.two_norm();

            const ctype denom = n*d;
            const ctype dist = n*(a-geo1.corner(f[0]));

            // if denominator is zero the segment in parallel to
            // the plane. If the distance is positive there is no intersection
            using std::abs;
            if (abs(denom) < eps)
            {
                if (dist > eps)
                    return false;
            }
            else // not parallel: compute line-plane intersection
            {
                const ctype t = -dist / denom;
                // if entering half space cut tfirst if t is larger
                using std::signbit;
                if (signbit(denom))
                {
                    if (t > tfirst)
                        tfirst = t;
                }
                // if exiting half space cut tlast if t is smaller
                else
                {
                    if (t < tlast)
                        tlast = t;
                }
                // there is no intersection of the interval is empty
                // use unscaled epsilon since t is in local coordinates
                if (tfirst > tlast - eps_)
                    return false;
            }
        }
        // If we made it until here an intersection exists. We also export
        // the intersections geometry now s(t) = a + t(b-a) in [tfirst, tlast]
        intersection = {std::vector<Point>({geo2.global(tfirst), geo2.global(tlast)})};
        return true;
    }
};

/*!
 * \ingroup Geometry
 * \brief A class for polyhedron--polygon intersection in 3d space
 */
template <class Geometry1, class Geometry2>
class GeometryIntersection<Geometry1, Geometry2, 3, 3, 2>
{
    enum { dimworld = 3 };
    enum { dim1 = 3 };
    enum { dim2 = 2 };

public:
    using ctype = typename Dune::PromotionTraits<typename Geometry1::ctype, typename Geometry2::ctype>::PromotedType;
    using Point = Dune::FieldVector<ctype, dimworld>;
    using IntersectionType = std::vector<std::array<Point, 3>>;

private:
    static constexpr ctype eps_ = 1.5e-7; // base epsilon for floating point comparisons
    using ReferenceElementsGeo1 = typename Dune::ReferenceElements<ctype, dim1>;
    using ReferenceElementsGeo2 = typename Dune::ReferenceElements<ctype, dim2>;

public:
    /*!
     * \brief Colliding segment and convex polyhedron
     * \note First we find the vertex candidates for the intersection region as follows:
     *       Add triangle vertices that are inside the tetrahedron
     *       Add tetrahedron vertices that are inside the triangle
     *       Add all intersection points of tetrahedron edges (codim 2) with the triangle (codim 0) (6*1 tests)
     *       Add all intersection points of triangle edges (codim 1) with tetrahedron faces (codim 1) (3*4 tests)
     *       Remove duplicate points from the list
     *       Compute the convex hull polygon
     *       Return a triangulation of that polygon as intersection
     * \param geo1/geo2 The geometries to intersect
     * \param intersection A triangulation of the intersection polygon
     */
    static bool intersection(const Geometry1& geo1, const Geometry2& geo2, IntersectionType& intersection)
    {
        static_assert(int(dimworld) == int(Geometry2::coorddimension),
                      "Can only collide geometries of same coordinate dimension");

        // the candidate intersection points
        std::vector<Point> points; points.reserve(10);

        // add 3d geometry corners that are inside the 2d geometry
        for (int i = 0; i < geo1.corners(); ++i)
            if (intersectsPointGeometry(geo1.corner(i), geo2))
                points.emplace_back(geo1.corner(i));

        // add 2d geometry corners that are inside the 3d geometry
        for (int i = 0; i < geo2.corners(); ++i)
            if (intersectsPointGeometry(geo2.corner(i), geo1))
                points.emplace_back(geo2.corner(i));

        // get some geometry types
        using PolyhedronFaceGeometry = Dune::MultiLinearGeometry<ctype, 2, dimworld>;
        using SegGeometry = Dune::MultiLinearGeometry<ctype, 1, dimworld>;

        const auto referenceElement1 = ReferenceElementsGeo1::general(geo1.type());
        const auto referenceElement2 = ReferenceElementsGeo2::general(geo2.type());

        // add intersection points of all polyhedron edges (codim dim-1) with the polygon
        for (int i = 0; i < referenceElement1.size(dim1-1); ++i)
        {
            const auto localEdgeGeom = referenceElement1.template geometry<dim1-1>(i);
            const auto p = geo1.global(localEdgeGeom.corner(0));
            const auto q = geo1.global(localEdgeGeom.corner(1));
            const auto segGeo = SegGeometry(Dune::GeometryTypes::line, std::vector<Point>{p, q});

            using PolySegTest = GeometryIntersection<Geometry2, SegGeometry>;
            typename PolySegTest::IntersectionType polySegIntersection;
            if (PolySegTest::template intersection<2>(geo2, segGeo, polySegIntersection))
                points.emplace_back(polySegIntersection[0]);
        }

        // add intersection points of all polygon faces (codim 1) with the polyhedron faces
        for (int i = 0; i < referenceElement1.size(1); ++i)
        {
            const auto faceGeo = [&]()
            {
                const auto localFaceGeo = referenceElement1.template geometry<1>(i);
                if (localFaceGeo.corners() == 4)
                {
                    const auto a = geo1.global(localFaceGeo.corner(0));
                    const auto b = geo1.global(localFaceGeo.corner(1));
                    const auto c = geo1.global(localFaceGeo.corner(2));
                    const auto d = geo1.global(localFaceGeo.corner(3));

                    return PolyhedronFaceGeometry(Dune::GeometryTypes::cube(2), std::vector<Point>{a, b, c, d});
                }
                else
                {
                    const auto a = geo1.global(localFaceGeo.corner(0));
                    const auto b = geo1.global(localFaceGeo.corner(1));
                    const auto c = geo1.global(localFaceGeo.corner(2));

                    return PolyhedronFaceGeometry(Dune::GeometryTypes::simplex(2), std::vector<Point>{a, b, c});
                }
            }();

            for (int j = 0; j < referenceElement2.size(1); ++j)
            {
                const auto localEdgeGeom = referenceElement2.template geometry<1>(j);
                const auto p = geo2.global(localEdgeGeom.corner(0));
                const auto q = geo2.global(localEdgeGeom.corner(1));

                const auto segGeo = SegGeometry(Dune::GeometryTypes::line, std::vector<Point>{p, q});

                using PolySegTest = GeometryIntersection<PolyhedronFaceGeometry, SegGeometry>;
                typename PolySegTest::IntersectionType polySegIntersection;
                if (PolySegTest::template intersection<2>(faceGeo, segGeo, polySegIntersection))
                    points.emplace_back(polySegIntersection[0]);
            }
        }

        // return if no intersection points were found
        if (points.empty()) return false;

        // remove duplicates
        const auto eps = (geo1.corner(0) - geo1.corner(1)).two_norm()*eps_;
        std::sort(points.begin(), points.end(), [&eps](const auto& a, const auto& b) -> bool
        {
            using std::abs;
            return (abs(a[0]-b[0]) > eps ? a[0] < b[0] : (abs(a[1]-b[1]) > eps ? a[1] < b[1] : (a[2] < b[2])));
        });

        auto removeIt = std::unique(points.begin(), points.end(), [&eps](const auto& a, const auto&b)
        {
            return (b-a).two_norm() < eps;
        });

        points.erase(removeIt, points.end());

        // return false if we don't have more than three unique points
        if (points.size() < 3) return false;

        // compute convex hull
        const auto convexHull = grahamConvexHull2d3d(points);
        assert(!convexHull.empty());

        // the intersections are the triangulation of the convex hull polygon
        intersection = triangulateConvexHull(convexHull);
        return true;
    }
};

/*!
 * \ingroup Geometry
 * \brief A class for polygon--segment intersection in 3d space
 */
template <class Geometry1, class Geometry2>
class GeometryIntersection<Geometry1, Geometry2, 3, 2, 1>
{
    enum { dimworld = 3 };
    enum { dim1 = 2 };
    enum { dim2 = 1 };

public:
    using ctype = typename Dune::PromotionTraits<typename Geometry1::ctype, typename Geometry2::ctype>::PromotedType;
    using Point = Dune::FieldVector<ctype, dimworld>;
    using IntersectionType = std::vector<Point>;

private:
    static constexpr ctype eps_ = 1.5e-7; // base epsilon for floating point comparisons
    using ReferenceElements = typename Dune::ReferenceElements<ctype, dim1>;

public:
    /*!
     *  \brief Colliding segment and convex polyhedron
     *  \note Algorithm based on the one from "Real-Time Collision Detection" by Christer Ericson,
     *        published by Morgan Kaufmann Publishers, (c) 2005 Elsevier Inc. (Chapter 5.3.6)
     * \param geo1/geo2 The geometries to intersect
     * \param is If the geometries collide, is holds the corner points of
     *        the intersection object in global coordinates.
     */
    template<int dimIntersection>
    static bool intersection(const Geometry1& geo1, const Geometry2& geo2, IntersectionType& is)
    {
        if (dimIntersection != 2)
            DUNE_THROW(Dune::NotImplemented, "Only simplex intersections are currently implemented!");

        static_assert(int(dimworld) == int(Geometry2::coorddimension), "Can only collide geometries of same coordinate dimension");

        const auto p = geo2.corner(0);
        const auto q = geo2.corner(1);

        const auto a = geo1.corner(0);
        const auto b = geo1.corner(1);
        const auto c = geo1.corner(2);

        if (geo1.corners() == 3)
            return intersection<dimIntersection>(a, b, c, p, q, is);

        else if (geo1.corners() == 4)
        {
            const auto d = geo1.corner(3);
            if (intersection<dimIntersection>(a, b, d, p, q, is)) return true;
            else if (intersection<dimIntersection>(a, d, c, p, q, is)) return true;
            else return false;
        }

        else
            DUNE_THROW(Dune::NotImplemented, "Collision of segment and geometry of type "
                          << geo1.type() << ", "<< geo1.corners() << " corners.");
    }

    // triangle--segment intersection with points as input
    template<int dimIntersection>
    static bool intersection(const Point& a, const Point& b, const Point& c,
                             const Point& p, const Point& q,
                             IntersectionType& is)
    {
        if (dimIntersection != 2)
            DUNE_THROW(Dune::NotImplemented, "Only simplex intersections are currently implemented!");

        const auto ab = b - a;
        const auto ac = c - a;
        const auto qp = p - q;

        // compute the triangle normal that defines the triangle plane
        const auto normal = crossProduct(ab, ac);

        // compute the denominator
        // if denom is 0 the segment is parallel and we can return
        const auto denom = normal*qp;
        const auto eps = eps_*ab.two_norm2()*qp.two_norm();
        using std::abs;
        if (abs(denom) < eps)
            return false;

        // compute intersection t value of pq with plane of triangle.
        // a segment intersects if and only if 0 <= t <= 1.
        const auto ap = p - a;
        const auto t = (ap*normal)/denom;
        if (t < 0.0 - eps_) return false;
        if (t > 1.0 + eps_) return false;

        // compute the barycentric coordinates and check if the intersection point
        // is inside the bounds of the triangle
        const auto e = crossProduct(qp, ap);
        const auto v = (ac*e)/denom;
        if (v < -eps_ || v > 1.0 + eps_) return false;
        const auto w = -(ab*e)/denom;
        if (w < -eps_ || v + w > 1.0 + eps_) return false;

        // Now we are sure there is an intersection points
        // Perform delayed division compute the last barycentric coordinate component
        const auto u = 1.0 - v - w;

        Point ip(0.0);
        ip.axpy(u, a);
        ip.axpy(v, b);
        ip.axpy(w, c);
        is = {ip};
        return true;
    }
};

} // end namespace Dumux

# endif
