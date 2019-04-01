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

#include <dune/common/deprecated.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/promotiontraits.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/math.hh>
#include <dumux/common/geometry/intersectspointgeometry.hh>
#include <dumux/common/geometry/grahamconvexhull.hh>

namespace Dumux {
namespace IntersectionPolicy {

//! Policy structure for point-like intersections
template<class ct, int dw>
struct PointPolicy
{
    using ctype = ct;
    using Point = Dune::FieldVector<ctype, dw>;

    static constexpr int dimWorld = dw;
    static constexpr int dimIntersection = 0;
    using Intersection = Point;
};

//! Policy structure for segment-like intersections
template<class ct, int dw>
struct SegmentPolicy
{
    using ctype = ct;
    using Point = Dune::FieldVector<ctype, dw>;
    using Segment = std::array<Point, 2>;

    static constexpr int dimWorld = dw;
    static constexpr int dimIntersection = 1;
    using Intersection = Segment;
};

//! Policy structure for polygon-like intersections
template<class ct, int dw>
struct PolygonPolicy
{
    using ctype = ct;
    using Point = Dune::FieldVector<ctype, dw>;
    using Polygon = std::vector<Point>;

    static constexpr int dimWorld = dw;
    static constexpr int dimIntersection = 2;
    using Intersection = Polygon;
};

//! default policy chooser class
template<class Geometry1, class Geometry2, int dim2 = Geometry2::mydimension>
class DefaultPolicyChooser
{
    static constexpr int dimworld = Geometry1::coorddimension;
    static_assert(dimworld == int(Geometry2::coorddimension), "Geometries must have the same coordinate dimension!");
    static_assert(dim2 == 1 || dim2 == 2, "No chooser class specialization available for given dimensions");
    using ctype = typename Dune::PromotionTraits<typename Geometry1::ctype, typename Geometry2::ctype>::PromotedType;
public:
    using type = std::conditional_t< dim2 == 2,
                                     PolygonPolicy<ctype, Geometry1::coorddimension>,
                                     SegmentPolicy<ctype, Geometry1::coorddimension> >;
};

//! Helper alias to define the default intersection policy
template<class Geometry1, class Geometry2>
using DefaultPolicy = typename DefaultPolicyChooser<Geometry1, Geometry2>::type;

} // end namespace IntersectionPolicy

/*!
 * \ingroup Geometry
 * \brief A class for geometry collision detection and intersection calculation
 * The class can be specialized for combinations of dimworld, dim1, dim2, where
 * dimworld is the world dimension embedding a grid of dim1 and a grid of dim2.
 */
template
<class Geometry1, class Geometry2,
 class Policy = IntersectionPolicy::DefaultPolicy<Geometry1, Geometry2>,
 int dimworld = Geometry1::coorddimension,
 int dim1 = Geometry1::mydimension,
 int dim2 = Geometry2::mydimension>
class GeometryIntersection
{
    static constexpr int dimWorld = Policy::dimWorld;

public:
    using ctype = typename Policy::ctype;
    using Point = typename Policy::Point;
    using Intersection = typename Policy::Intersection;

    //! Deprecated alias, will be removed after 3.1
    using IntersectionType DUNE_DEPRECATED_MSG("Please use Intersection instead") = Intersection;

    //! Determine if the two geometries intersect and compute the intersection geometry
    static bool intersection(const Geometry1& geo1, const Geometry2& geo2, Intersection& intersection)
    {
        static_assert(int(dim1) >= int(dim2), "Geometries must be ordered such that dim1 >= dim2");
        static_assert(dimworld == Geometry2::coorddimension, "Can only intersect geometries of same coordinate dimension");
        DUNE_THROW(Dune::NotImplemented, "Geometry intersection detection with intersection computation not implemented for dimworld = "
                                          << dimworld << ", dim1 = " << dim1 << ", dim2 = " << dim2);
    }
};

/*!
 * \ingroup Geometry
 * \brief A class for polygon--segment intersection in 2d space
 */
template <class Geometry1, class Geometry2, class Policy>
class GeometryIntersection<Geometry1, Geometry2, Policy, 2, 2, 1>
{
    enum { dimworld = 2 };
    enum { dim1 = 2 };
    enum { dim2 = 1 };

public:
    using ctype = typename Policy::ctype;
    using Point = typename Policy::Point;
    using Intersection = typename Policy::Intersection;

    //! Deprecated alias, will be removed after 3.1
    using IntersectionType DUNE_DEPRECATED_MSG("Please use Intersection instead") = Intersection;

private:
    static constexpr ctype eps_ = 1.5e-7; // base epsilon for floating point comparisons
    using ReferenceElements = typename Dune::ReferenceElements<ctype, dim1>;

public:
    /*!
     * \brief Colliding segment and convex polygon
     * \param geo1/geo2 The geometries to intersect
     * \param intersection If the geometries collide intersection holds the
     *        corner points of the intersection object in global coordinates.
     * \note This overload is used when segment-like intersections are seeked.
     */
    template<class P = Policy, std::enable_if_t<P::dimIntersection == 1, int> = 0>
    static bool intersection(const Geometry1& geo1, const Geometry2& geo2, Intersection& intersection)
    {
        static_assert(int(dimworld) == int(Geometry2::coorddimension), "Can only collide geometries of same coordinate dimension");

        ctype tfirst, tlast;
        if (intersect_(geo1, geo2, tfirst, tlast))
        {
            // the intersection exists. Export the intersections geometry now:
            // s(t) = a + t(b-a) in [tfirst, tlast]
            intersection = Intersection({geo2.global(tfirst), geo2.global(tlast)});
            return true;
        }

        return false;
    }

    /*!
     * \brief Colliding segment and convex polygon
     * \param geo1/geo2 The geometries to intersect
     * \param intersection If the geometries collide intersection holds the
     *        corner points of the intersection object in global coordinates.
     * \note This overload is used when point-like intersections are seeked.
     */
    template<class P = Policy, std::enable_if_t<P::dimIntersection == 0, int> = 0>
    static bool intersection(const Geometry1& geo1, const Geometry2& geo2, typename P::Intersection& intersection)
    {
        static_assert(int(dimworld) == int(Geometry2::coorddimension), "Can only collide geometries of same coordinate dimension");

        ctype tfirst, tlast;
        if (!intersect_(geo1, geo2, tfirst, tlast))
        {
            // if start & end point are same, the intersection is a point
            using std::abs;
            if ( tfirst > tlast - eps_ )
            {
                intersection = Intersection(geo2.global(tfirst));
                return true;
            }
        }

        return false;
    }

private:
    /*!
     * \brief Obtain local coordinates of start/end point of the intersecting segment.
     */
     static bool intersect_(const Geometry1& geo1, const Geometry2& geo2, ctype& tfirst, ctype& tlast)
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
         tfirst = 0.0;
         tlast = 1.0;

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

         // If we made it until here an intersection exists.
         return true;
     }
};

/*!
 * \ingroup Geometry
 * \brief A class for polyhedron--segment intersection in 3d space
 */
template <class Geometry1, class Geometry2, class Policy>
class GeometryIntersection<Geometry1, Geometry2, Policy, 3, 3, 1>
{
    enum { dimworld = 3 };
    enum { dim1 = 3 };
    enum { dim2 = 1 };

public:
    using ctype = typename Policy::ctype;
    using Point = typename Policy::Point;
    using Intersection = typename Policy::Intersection;

    //! Deprecated alias, will be removed after 3.1
    using IntersectionType DUNE_DEPRECATED_MSG("Please use Intersection instead") = Intersection;

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
     * \note This overload is used when segment-like intersections are seeked.
     */
    template<class P = Policy, std::enable_if_t<P::dimIntersection == 1, int> = 0>
    static bool intersection(const Geometry1& geo1, const Geometry2& geo2, Intersection& intersection)
    {
        static_assert(int(dimworld) == int(Geometry2::coorddimension), "Can only collide geometries of same coordinate dimension");

        ctype tfirst, tlast;
        if (intersect_(geo1, geo2, tfirst, tlast))
        {
            // Intersection exists. Export the intersections geometry now:
            // s(t) = a + t(b-a) in [tfirst, tlast]
            intersection = Intersection({geo2.global(tfirst), geo2.global(tlast)});
            return true;
        }

        return false;
    }

    /*!
     *  \brief Colliding segment and convex polyhedron
     *  \note Algorithm based on the one from "Real-Time Collision Detection" by Christer Ericson,
     *        published by Morgan Kaufmann Publishers, (c) 2005 Elsevier Inc.
     *        Basis is the theorem that for any two non-intersecting convex polyhedrons
     *        a separating plane exists.
     * \param geo1/geo2 The geometries to intersect
     * \param intersection If the geometries collide intersection holds the corner points of
     *        the intersection object in global coordinates.
     * \note This overload is used when point-like intersections are seeked.
     */
    template<class P = Policy, std::enable_if_t<P::dimIntersection == 0, int> = 0>
    static bool intersection(const Geometry1& geo1, const Geometry2& geo2, Intersection& intersection)
    {
        static_assert(int(dimworld) == int(Geometry2::coorddimension), "Can only collide geometries of same coordinate dimension");

        ctype tfirst, tlast;
        if (intersect_(geo1, geo2, tfirst, tlast))
        {
            // if start & end point are same, the intersection is a point
            using std::abs;
            if ( abs(tfirst - tlast) < eps_ )
            {
                intersection = Intersection(geo2.global(tfirst));
                return true;
            }
        }

        return false;
    }

private:
    /*!
     * \brief Obtain local coordinates of start/end point of the intersecting segment.
     */
     static bool intersect_(const Geometry1& geo1, const Geometry2& geo2, ctype& tfirst, ctype& tlast)
     {
         const auto a = geo2.corner(0);
         const auto b = geo2.corner(1);
         const auto d = b - a;

         // The initial interval is the whole segment
         // afterward we start clipping the interval
         // by the planes decribed by the facet
         tfirst = 0.0;
         tlast = 1.0;

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

         // If we made it until here an intersection exists.
         return true;
     }
};

/*!
 * \ingroup Geometry
 * \brief A class for polyhedron--polygon intersection in 3d space
 */
template <class Geometry1, class Geometry2, class Policy>
class GeometryIntersection<Geometry1, Geometry2, Policy, 3, 3, 2>
{
    enum { dimworld = 3 };
    enum { dim1 = 3 };
    enum { dim2 = 2 };

public:
    using ctype = typename Policy::ctype;
    using Point = typename Policy::Point;
    using Intersection = typename Policy::Intersection;

    //! Deprecated alias, will be removed after 3.1
    using IntersectionType DUNE_DEPRECATED_MSG("Please use Intersection instead") = Intersection;

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
     * \note This overload is used when polygon like intersections are seeked
     */
    template<class P = Policy, std::enable_if_t<P::dimIntersection == 2, int> = 0>
    static bool intersection(const Geometry1& geo1, const Geometry2& geo2, Intersection& intersection)
    {
        static_assert(int(dimworld) == int(Geometry2::coorddimension), "Can only collide geometries of same coordinate dimension");

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

        // intersection policy for point-like intersections (used later)
        using PointPolicy = IntersectionPolicy::PointPolicy<ctype, dimworld>;

        // add intersection points of all polyhedron edges (codim dim-1) with the polygon
        for (int i = 0; i < referenceElement1.size(dim1-1); ++i)
        {
            const auto localEdgeGeom = referenceElement1.template geometry<dim1-1>(i);
            const auto p = geo1.global(localEdgeGeom.corner(0));
            const auto q = geo1.global(localEdgeGeom.corner(1));
            const auto segGeo = SegGeometry(Dune::GeometryTypes::line, std::vector<Point>{p, q});

            using PolySegTest = GeometryIntersection<Geometry2, SegGeometry, PointPolicy>;
            typename PolySegTest::Intersection polySegIntersection;
            if (PolySegTest::intersection(geo2, segGeo, polySegIntersection))
                points.emplace_back(polySegIntersection);
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

                using PolySegTest = GeometryIntersection<PolyhedronFaceGeometry, SegGeometry, PointPolicy>;
                typename PolySegTest::Intersection polySegIntersection;
                if (PolySegTest::intersection(faceGeo, segGeo, polySegIntersection))
                    points.emplace_back(polySegIntersection);
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

        // intersection polygon is convex hull of above points
        intersection = grahamConvexHull2d3d(points);
        assert(!intersection.empty());

        return true;
    }

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
     * \todo implement overloads for segment or point-like intersections
     */
    template<class P = Policy, std::enable_if_t<P::dimIntersection != 2, int> = 0>
    static bool intersection(const Geometry1& geo1, const Geometry2& geo2, Intersection& intersection)
    {
        static_assert(int(dimworld) == int(Geometry2::coorddimension), "Can only collide geometries of same coordinate dimension");
        DUNE_THROW(Dune::NotImplemented, "Polyhedron-polygon intersection detection only implemented for polygon-like intersections");
    }
};

/*!
 * \ingroup Geometry
 * \brief A class for polygon--segment intersection in 3d space
 */
template <class Geometry1, class Geometry2, class Policy>
class GeometryIntersection<Geometry1, Geometry2, Policy, 3, 2, 1>
{
    enum { dimworld = 3 };
    enum { dim1 = 2 };
    enum { dim2 = 1 };

public:
    using ctype = typename Policy::ctype;
    using Point = typename Policy::Point;
    using Intersection = typename Policy::Intersection;

    //! Deprecated alias, will be removed after 3.1
    using IntersectionType DUNE_DEPRECATED_MSG("Please use Intersection instead") = Intersection;

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
     * \note This overload is used when point-like intersections are seeked
     */
    template<class P = Policy, std::enable_if_t<P::dimIntersection == 0, int> = 0>
    static bool intersection(const Geometry1& geo1, const Geometry2& geo2, Intersection& is)
    {
        static_assert(int(dimworld) == int(Geometry2::coorddimension), "Can only collide geometries of same coordinate dimension");

        const auto p = geo2.corner(0);
        const auto q = geo2.corner(1);

        const auto a = geo1.corner(0);
        const auto b = geo1.corner(1);
        const auto c = geo1.corner(2);

        if (geo1.corners() == 3)
            return intersect_<Policy>(a, b, c, p, q, is);

        else if (geo1.corners() == 4)
        {
            Intersection is1, is2;
            bool hasSegment1, hasSegment2;

            const auto d = geo1.corner(3);
            const bool intersects1 = intersect_<Policy>(a, b, d, p, q, is1, hasSegment1);
            const bool intersects2 = intersect_<Policy>(a, d, c, p, q, is2, hasSegment2);

            if (hasSegment1 || hasSegment2)
                return false;

            if (intersects1) { is = is1; return true; }
            if (intersects2) { is = is2; return true; }

            return false;
        }

        else
            DUNE_THROW(Dune::NotImplemented, "Collision of segment and geometry of type "
                          << geo1.type() << ", "<< geo1.corners() << " corners.");
    }

    /*!
     *  \brief Colliding segment and convex polyhedron
     *  \note Algorithm based on the one from "Real-Time Collision Detection" by Christer Ericson,
     *        published by Morgan Kaufmann Publishers, (c) 2005 Elsevier Inc. (Chapter 5.3.6)
     * \param geo1/geo2 The geometries to intersect
     * \param is If the geometries collide, is holds the corner points of
     *        the intersection object in global coordinates.
     * \note This overload is used when point-like intersections are seeked
     */
    template<class P = Policy, std::enable_if_t<P::dimIntersection == 0, int> = 0>
    DUNE_DEPRECATED_MSG("Please use intersection(triangle, segment, ...) instead")
    static bool intersection(const Point& a, const Point& b, const Point& c,
                             const Point& p, const Point& q,
                             Intersection& is)
    {
        return intersect_<Policy>(a, b, c, p, q, is);
    }

private:
    // triangle--segment intersection with points as input
    template<class P = Policy, std::enable_if_t<P::dimIntersection == 0, int> = 0>
    static bool intersect_(const Point& a, const Point& b, const Point& c,
                           const Point& p, const Point& q,
                           Intersection& is)
    {
        bool hasSegment;
        return intersect_(a, b, c, p, q, is, hasSegment);
    }

    /*!
     * \brief triangle--segment point-like intersection with points as input. Also
     *        stores if a segment-like intersection was found in the provided boolean.
     */
    template<class P = Policy, std::enable_if_t<P::dimIntersection == 0, int> = 0>
    static bool intersect_(const Point& a, const Point& b, const Point& c,
                           const Point& p, const Point& q,
                           Intersection& is, bool& hasSegmentIs)
    {
        hasSegmentIs = false;

        const auto ab = b - a;
        const auto ac = c - a;
        const auto qp = p - q;

        // compute the triangle normal that defines the triangle plane
        const auto normal = crossProduct(ab, ac);

        // compute the denominator
        // if denom is 0 the segment is parallel and we can return
        const auto denom = normal*qp;
        const auto abnorm2 = ab.two_norm2();
        const auto eps = eps_*abnorm2*qp.two_norm();
        using std::abs;
        if (abs(denom) < eps)
        {
            const auto pa = a - p;
            const auto denom2 = normal*pa;
            if (abs(denom2) > eps_*pa.two_norm()*abnorm2)
                return false;

            // if we get here, we are in-plane. Check if a
            // segment-like intersection with geometry 1 exists.
            using SegmentPolicy = typename IntersectionPolicy::SegmentPolicy<ctype, dimworld>;
            using Triangle = Dune::AffineGeometry<ctype, 2, dimworld>;
            using Segment = Dune::AffineGeometry<ctype, 1, dimworld>;
            using SegmentIntersectionAlgorithm = GeometryIntersection<Triangle, Segment, SegmentPolicy>;
            using SegmentIntersectionType = typename SegmentIntersectionAlgorithm::Intersection;
            SegmentIntersectionType segmentIs;

            Triangle triangle(Dune::GeometryTypes::simplex(2), std::array<Point, 3>({a, b, c}));
            Segment segment(Dune::GeometryTypes::simplex(1), std::array<Point, 2>({p, q}));
            if (SegmentIntersectionAlgorithm::intersection(triangle, segment, segmentIs))
            {
                hasSegmentIs = true;
                return false;
            }

            // We are in-plane. A point-like
            // intersection can only be on the edges
            for (const auto& ip : {p, q})
            {
                if ( intersectsPointSimplex(ip, a, b)
                     || intersectsPointSimplex(ip, b, c)
                     || intersectsPointSimplex(ip, c, a) )
                {
                    is = ip;
                    return true;
                }
            }

            return false;
        }

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
        is = ip;
        return true;
    }
};

} // end namespace Dumux

# endif
