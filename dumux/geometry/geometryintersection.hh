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

#include <tuple>

#include <dune/common/exceptions.hh>
#include <dune/common/promotiontraits.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/math.hh>
#include <dumux/geometry/intersectspointgeometry.hh>
#include <dumux/geometry/grahamconvexhull.hh>
#include <dumux/geometry/boundingboxtree.hh>

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

//! Policy structure for polyhedron-like intersections
template<class ct, int dw>
struct PolyhedronPolicy
{
    using ctype = ct;
    using Point = Dune::FieldVector<ctype, dw>;
    using Polyhedron = std::vector<Point>;

    static constexpr int dimWorld = dw;
    static constexpr int dimIntersection = 3;
    using Intersection = Polyhedron;
};

//! default policy chooser class
template<class Geometry1, class Geometry2>
class DefaultPolicyChooser
{
    static constexpr int dimworld = Geometry1::coorddimension;
    static constexpr int isDim = std::min( int(Geometry1::mydimension), int(Geometry2::mydimension) );
    static_assert(dimworld == int(Geometry2::coorddimension),
                  "Geometries must have the same coordinate dimension!");
    static_assert(int(Geometry1::mydimension) <= 3 && int(Geometry2::mydimension) <= 3,
                  "Geometries must have dimension 3 or less.");
    using ctype = typename Dune::PromotionTraits<typename Geometry1::ctype, typename Geometry2::ctype>::PromotedType;

    using DefaultPolicies = std::tuple<PointPolicy<ctype, dimworld>,
                                       SegmentPolicy<ctype, dimworld>,
                                       PolygonPolicy<ctype, dimworld>,
                                       PolyhedronPolicy<ctype, dimworld>>;
public:
    using type = std::tuple_element_t<isDim, DefaultPolicies>;
};

//! Helper alias to define the default intersection policy
template<class Geometry1, class Geometry2>
using DefaultPolicy = typename DefaultPolicyChooser<Geometry1, Geometry2>::type;

} // end namespace IntersectionPolicy

namespace Detail {

/*!
 * \ingroup Geometry
 * \brief Algorithm to find segment-like intersections of a polgon/polyhedron with a
 *        segment. The result is stored in the form of the local coordinates tfirst
 *        and tlast on the segment geo1.
 * \param geo1 the first geometry
 * \param geo2 the second geometry (dimGeo2 < dimGeo1)
 * \param baseEps the base epsilon used for floating point comparisons
 * \param tfirst stores the local coordinate of beginning of intersection segment
 * \param tlast stores the local coordinate of the end of intersection segment
 * \param getFacetCornerIndices Function to obtain the facet corner indices on geo1
 * \param computeNormal Function to obtain the normal vector on a facet on geo1
 */
template< class Geo1, class Geo2, class ctype,
          class GetFacetCornerIndices, class ComputeNormalFunction >
bool computeSegmentIntersection(const Geo1& geo1, const Geo2& geo2, ctype baseEps,
                                ctype& tfirst, ctype& tlast,
                                const GetFacetCornerIndices& getFacetCornerIndices,
                                const ComputeNormalFunction& computeNormal)
{
    static_assert(int(Geo2::mydimension) == 1, "Geometry2 must be a segment");
    static_assert(int(Geo1::mydimension) > int(Geo2::mydimension),
                  "Dimension of Geometry1 must be higher than that of Geometry2");

    const auto a = geo2.corner(0);
    const auto b = geo2.corner(1);
    const auto d = b - a;

    // The initial interval is the whole segment.
    // Afterwards we start clipping the interval
    // at the intersections with the facets of geo1
    tfirst = 0.0;
    tlast = 1.0;

    const auto facets = getFacetCornerIndices(geo1);
    for (const auto& f : facets)
    {
        const auto n = computeNormal(f);

        const auto c0 = geo1.corner(f[0]);
        const ctype denom = n*d;
        const ctype dist = n*(a-c0);

        // use first edge of the facet to scale eps
        const auto edge1 = geo1.corner(f[1]) - geo1.corner(f[0]);
        const ctype eps = baseEps*edge1.two_norm();

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
            if (tfirst > tlast - baseEps)
                return false;
        }
    }

    // If we made it until here an intersection exists.
    return true;
}

} // end namespace Detail

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

    //! Determine if the two geometries intersect and compute the intersection geometry
    static bool intersection(const Geometry1& geo1, const Geometry2& geo2, Intersection& intersection)
    {
        static_assert(dimworld == Geometry2::coorddimension, "Can only intersect geometries of same coordinate dimension");
        DUNE_THROW(Dune::NotImplemented, "Geometry intersection detection with intersection computation not implemented for dimworld = "
                                          << dimworld << ", dim1 = " << dim1 << ", dim2 = " << dim2);
    }
};

/*!
 * \ingroup Geometry
 * \brief A class for segment--segment intersection in 2d space
 */
template <class Geometry1, class Geometry2, class Policy>
class GeometryIntersection<Geometry1, Geometry2, Policy, 2, 1, 1>
{
    enum { dimworld = 2 };
    enum { dim1 = 1 };
    enum { dim2 = 1 };

    // base epsilon for floating point comparisons
    static constexpr typename Policy::ctype eps_ = 1.5e-7;

public:
    using ctype = typename Policy::ctype;
    using Point = typename Policy::Point;
    using Intersection = typename Policy::Intersection;

    /*!
     * \brief Colliding two segments
     * \param geo1/geo2 The geometries to intersect
     * \param intersection The intersection point
     * \note This overload is used when point-like intersections are seeked
     */
    template<class P = Policy, std::enable_if_t<P::dimIntersection == 0, int> = 0>
    static bool intersection(const Geometry1& geo1, const Geometry2& geo2, Intersection& intersection)
    {
        static_assert(int(dimworld) == int(Geometry2::coorddimension), "Can only collide geometries of same coordinate dimension");

        const auto v1 = geo1.corner(1) - geo1.corner(0);
        const auto v2 = geo2.corner(1) - geo2.corner(0);
        const auto ac = geo2.corner(0) - geo1.corner(0);

        auto n2 = Point({-1.0*v2[1], v2[0]});
        n2 /= n2.two_norm();

        // compute distance of first corner on geo2 to line1
        const auto dist12 = n2*ac;

        // first check parallel segments
        using std::abs;
        using std::sqrt;

        const auto v1Norm2 = v1.two_norm2();
        const auto eps = eps_*sqrt(v1Norm2);
        const auto eps2 = eps_*v1Norm2;

        const auto sp = n2*v1;
        if (abs(sp) < eps)
        {
            if (abs(dist12) > eps)
                return false;

            // point intersection can only be one of the corners
            if ( (v1*v2) < 0.0 )
            {
                if ( ac.two_norm2() < eps2 )
                { intersection = geo2.corner(0); return true; }

                if ( (geo2.corner(1) - geo1.corner(1)).two_norm2() < eps2 )
                { intersection = geo2.corner(1); return true; }
            }
            else
            {
                if ( (geo2.corner(1) - geo1.corner(0)).two_norm2() < eps2 )
                { intersection = geo2.corner(1); return true; }

                if ( (geo2.corner(0) - geo1.corner(1)).two_norm2() < eps2 )
                { intersection = geo2.corner(0); return true; }
            }

            // no intersection
            return false;
        }

        // intersection point on v1 in local coords
        const auto t1 = dist12 / sp;

        // check if the local coords are valid
        if (t1 < -1.0*eps_ || t1 > 1.0 + eps_)
            return false;

        // compute global coordinates
        auto isPoint = geo1.global(t1);

        // check if point is in bounding box of v2
        const auto c = geo2.corner(0);
        const auto d = geo2.corner(1);

        using std::min; using std::max;
        std::array<ctype, 4> bBox({ min(c[0], d[0]), min(c[1], d[1]), max(c[0], d[0]), max(c[1], d[1]) });

        if ( intersectsPointBoundingBox(isPoint, bBox.data()) )
        {
            intersection = std::move(isPoint);
            return true;
        }

        return false;
    }

    /*!
     * \brief Colliding two segments
     * \param geo1/geo2 The geometries to intersect
     * \param intersection Container to store corners of intersection segment
     * \note this overload is used when segment-like intersections are seeked
     */
    template<class P = Policy, std::enable_if_t<P::dimIntersection == 1, int> = 0>
    static bool intersection(const Geometry1& geo1, const Geometry2& geo2, Intersection& intersection)
    {
        static_assert(int(dimworld) == int(Geometry2::coorddimension), "Can only collide geometries of same coordinate dimension");

        const auto& a = geo1.corner(0);
        const auto& b = geo1.corner(1);
        const auto ab = b - a;

        const auto& p = geo2.corner(0);
        const auto& q = geo2.corner(1);
        const auto pq = q - p;
        Dune::FieldVector<ctype, dimworld> n2{-pq[1], pq[0]};

        using std::max;
        const auto abNorm2 = ab.two_norm2();
        const auto pqNorm2 = pq.two_norm2();
        const auto eps2 = eps_*max(abNorm2, pqNorm2);

        // non-parallel segments do not intersect in a segment.
        using std::abs;
        if (abs(n2*ab) > eps2)
            return false;

        // check if the segments lie on the same line
        const auto ap = p - a;
        if (abs(ap*n2) > eps2)
            return false;

        // compute scaled local coordinates of corner 1/2 of segment2 on segment1
        auto t1 = ab*ap;
        auto t2 = ab*(q - a);

        using std::swap;
        if (t1 > t2)
            swap(t1, t2);

        using std::clamp;
        t1 = clamp(t1, 0.0, abNorm2);
        t2 = clamp(t2, 0.0, abNorm2);

        if (abs(t2-t1) < eps2)
            return false;

        intersection = Intersection({geo1.global(t1/abNorm2), geo1.global(t2/abNorm2)});
        return true;
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

private:
    static constexpr ctype eps_ = 1.5e-7; // base epsilon for floating point comparisons

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
         // lambda to obtain the facet corners on geo1
         auto getFacetCorners = [] (const Geometry1& geo1)
         {
             std::vector< std::array<int, 2> > facetCorners;
             switch (geo1.corners())
             {
                 case 4: // quadrilateral
                     facetCorners = {{0, 2}, {3, 1}, {1, 0}, {2, 3}};
                     break;
                 case 3: // triangle
                     facetCorners = {{1, 0}, {0, 2}, {2, 1}};
                     break;
                 default:
                     DUNE_THROW(Dune::NotImplemented, "Collision of segment and geometry of type "
                                    << geo1.type() << " with "<< geo1.corners() << " corners.");
             }

             return facetCorners;
         };

         // lambda to obtain the normal on a facet
         const auto center1 = geo1.center();
         auto computeNormal = [&geo1, &center1] (const std::array<int, 2>& facetCorners)
         {
             const auto c0 = geo1.corner(facetCorners[0]);
             const auto c1 = geo1.corner(facetCorners[1]);
             const auto edge = c1 - c0;

             Dune::FieldVector<ctype, dimworld> n({-1.0*edge[1], edge[0]});
             n /= n.two_norm();

             // orientation of the element is not known
             // make sure normal points outwards of element
             if ( n*(center1-c0) > 0.0 )
                 n *= -1.0;

             return n;
         };

         return Detail::computeSegmentIntersection(geo1, geo2, eps_, tfirst, tlast, getFacetCorners, computeNormal);
     }
};

/*!
 * \ingroup Geometry
 * \brief A class for segment--polygon intersection in 2d space
 */
template <class Geometry1, class Geometry2, class Policy>
class GeometryIntersection<Geometry1, Geometry2, Policy, 2, 1, 2>
: public GeometryIntersection<Geometry2, Geometry1, Policy, 2, 2, 1>
{
    using Base = GeometryIntersection<Geometry2, Geometry1, Policy, 2, 2, 1>;

public:
    /*!
     * \brief Colliding segment and convex polygon
     * \param geo1/geo2 The geometries to intersect
     * \param intersection If the geometries collide intersection holds the
     *        corner points of the intersection object in global coordinates.
     * \note This forwards to the polygon-segment specialization with swapped arguments.
     */
    template<class P = Policy>
    static bool intersection(const Geometry1& geo1, const Geometry2& geo2, typename Base::Intersection& intersection)
    {
        return Base::intersection(geo2, geo1, intersection);
    }
};

/*!
 * \ingroup Geometry
 * \brief A class for polygon--polygon intersection in 2d space
 */
template <class Geometry1, class Geometry2, class Policy>
class GeometryIntersection<Geometry1, Geometry2, Policy, 2, 2, 2>
{
    enum { dimworld = 2 };
    enum { dim1 = 2 };
    enum { dim2 = 2 };

public:
    using ctype = typename Policy::ctype;
    using Point = typename Policy::Point;
    using Intersection = typename Policy::Intersection;

private:
    static constexpr ctype eps_ = 1.5e-7; // base epsilon for floating point comparisons

public:
    /*!
     * \brief Colliding two polygons
     * \note First we find the vertex candidates for the intersection region as follows:
     *       Add polygon vertices that are inside the other polygon
     *       Add intersections of polygon edges
     *       Remove duplicate points from the list
     *       Compute the convex hull polygon
     *       Return a triangulation of that polygon as intersection
     * \param geo1/geo2 The geometries to intersect
     * \param intersection Container to store the corner points of the polygon (as convex hull)
     * \note This overload is used when polygon like intersections are seeked
     */
    template<class P = Policy, std::enable_if_t<P::dimIntersection == 2, int> = 0>
    static bool intersection(const Geometry1& geo1, const Geometry2& geo2, Intersection& intersection)
    {
        static_assert(int(dimworld) == int(Geometry2::coorddimension), "Can only collide geometries of same coordinate dimension");

        // the candidate intersection points
        std::vector<Point> points; points.reserve(6);

        // add polygon1 corners that are inside polygon2
        for (int i = 0; i < geo1.corners(); ++i)
            if (intersectsPointGeometry(geo1.corner(i), geo2))
                points.emplace_back(geo1.corner(i));

        const auto numPoints1 = points.size();
        if (numPoints1 != geo1.corners())
        {
            // add polygon2 corners that are inside polygon1
            for (int i = 0; i < geo2.corners(); ++i)
                if (intersectsPointGeometry(geo2.corner(i), geo1))
                    points.emplace_back(geo2.corner(i));

            if (points.empty())
                return false;

            if (points.size() - numPoints1 != geo2.corners())
            {
                const auto refElement1 = referenceElement(geo1);
                const auto refElement2 = referenceElement(geo2);

                // add intersections of edges
                using SegGeometry = Dune::MultiLinearGeometry<ctype, 1, dimworld>;
                using PointPolicy = IntersectionPolicy::PointPolicy<ctype, dimworld>;
                for (int i = 0; i < refElement1.size(dim1-1); ++i)
                {
                    const auto localEdgeGeom1 = refElement1.template geometry<dim1-1>(i);
                    const auto edge1 = SegGeometry( Dune::GeometryTypes::line,
                                                    std::vector<Point>( {geo1.global(localEdgeGeom1.corner(0)),
                                                                         geo1.global(localEdgeGeom1.corner(1))} ));

                    for (int j = 0; j < refElement2.size(dim2-1); ++j)
                    {
                        const auto localEdgeGeom2 = refElement2.template geometry<dim2-1>(j);
                        const auto edge2 = SegGeometry( Dune::GeometryTypes::line,
                                                        std::vector<Point>( {geo2.global(localEdgeGeom2.corner(0)),
                                                                             geo2.global(localEdgeGeom2.corner(1))} ));

                        using EdgeTest = GeometryIntersection<SegGeometry, SegGeometry, PointPolicy>;
                        typename EdgeTest::Intersection edgeIntersection;
                        if (EdgeTest::intersection(edge1, edge2, edgeIntersection))
                            points.emplace_back(edgeIntersection);
                    }
                }
            }
        }

        if (points.empty())
            return false;

        // remove duplicates
        const auto eps = (geo1.corner(0) - geo1.corner(1)).two_norm()*eps_;
        std::sort(points.begin(), points.end(), [&eps](const auto& a, const auto& b) -> bool
        {
            using std::abs;
            return (abs(a[0]-b[0]) > eps ? a[0] < b[0] : a[1] < b[1]);
        });

        auto removeIt = std::unique(points.begin(), points.end(), [&eps](const auto& a, const auto&b)
        {
            return (b-a).two_norm() < eps;
        });

        points.erase(removeIt, points.end());

        // return false if we don't have at least three unique points
        if (points.size() < 3)
            return false;

        // intersection polygon is convex hull of above points
        intersection = grahamConvexHull<2>(points);
        assert(!intersection.empty());
        return true;
    }

    /*!
     * \brief Colliding two polygons
     * \param geo1/geo2 The geometries to intersect
     * \param intersection Container to store the corners of intersection segment
     * \note this overload is used when segment-like intersections are seeked
     * \todo implement this query
     */
    template<class P = Policy, std::enable_if_t<P::dimIntersection == 1, int> = 0>
    static bool intersection(const Geometry1& geo1, const Geometry2& geo2, Intersection& intersection)
    {
        static_assert(int(dimworld) == int(Geometry2::coorddimension), "Can only collide geometries of same coordinate dimension");
        DUNE_THROW(Dune::NotImplemented, "Polygon-polygon intersection detection for segment-like intersections");
    }

    /*!
     * \brief Colliding two polygons
     * \param geo1/geo2 The geometries to intersect
     * \param intersection The intersection point
     * \note this overload is used when point-like intersections are seeked
     * \todo implement this query
     */
    template<class P = Policy, std::enable_if_t<P::dimIntersection == 0, int> = 0>
    static bool intersection(const Geometry1& geo1, const Geometry2& geo2, Intersection& intersection)
    {
        static_assert(int(dimworld) == int(Geometry2::coorddimension), "Can only collide geometries of same coordinate dimension");
        DUNE_THROW(Dune::NotImplemented, "Polygon-polygon intersection detection for touching points");
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

private:
    static constexpr ctype eps_ = 1.5e-7; // base epsilon for floating point comparisons

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
         // lambda to obtain facet corners on geo1
         auto getFacetCorners = [] (const Geometry1& geo1)
         {
             std::vector< std::vector<int> > facetCorners;
             // sort facet corner so that normal n = (p1-p0)x(p2-p0) always points outwards
             switch (geo1.corners())
             {
                 // todo compute intersection geometries
                 case 8: // hexahedron
                     facetCorners = {{2, 0, 6, 4}, {1, 3, 5, 7}, {0, 1, 4, 5},
                               {3, 2, 7, 6}, {1, 0, 3, 2}, {4, 5, 6, 7}};
                     break;
                 case 4: // tetrahedron
                     facetCorners = {{1, 0, 2}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3}};
                     break;
                 default:
                     DUNE_THROW(Dune::NotImplemented, "Collision of segment and geometry of type "
                                    << geo1.type() << ", "<< geo1.corners() << " corners.");
             }

             return facetCorners;
         };

         // lambda to obtain the normal on a facet
         auto computeNormal = [&geo1] (const std::vector<int>& facetCorners)
         {
             const auto v0 = geo1.corner(facetCorners[1]) - geo1.corner(facetCorners[0]);
             const auto v1 = geo1.corner(facetCorners[2]) - geo1.corner(facetCorners[0]);

             auto n = crossProduct(v0, v1);
             n /= n.two_norm();

             return n;
         };

         return Detail::computeSegmentIntersection(geo1, geo2, eps_, tfirst, tlast, getFacetCorners, computeNormal);
     }
};

/*!
 * \ingroup Geometry
 * \brief A class for segment--polyhedron intersection in 3d space
 */
template <class Geometry1, class Geometry2, class Policy>
class GeometryIntersection<Geometry1, Geometry2, Policy, 3, 1, 3>
: public GeometryIntersection<Geometry2, Geometry1, Policy, 3, 3, 1>
{
    using Base = GeometryIntersection<Geometry2, Geometry1, Policy, 3, 3, 1>;
public:
    /*!
     * \brief Colliding segment and convex polyhedron
     * \param geo1/geo2 The geometries to intersect
     * \param intersection If the geometries collide intersection holds the
     *        corner points of the intersection object in global coordinates.
     * \note This forwards to the polyhedron-segment specialization with swapped arguments.
     */
    template<class P = Policy>
    static bool intersection(const Geometry1& geo1, const Geometry2& geo2, typename Base::Intersection& intersection)
    {
        return Base::intersection(geo2, geo1, intersection);
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

private:
    static constexpr ctype eps_ = 1.5e-7; // base epsilon for floating point comparisons

public:
    /*!
     * \brief Colliding polygon and convex polyhedron
     * \note First we find the vertex candidates for the intersection region as follows:
     *       Add triangle vertices that are inside the tetrahedron
     *       Add tetrahedron vertices that are inside the triangle
     *       Add all intersection points of tetrahedron edges (codim 2) with the triangle (codim 0) (6*1 tests)
     *       Add all intersection points of triangle edges (codim 1) with tetrahedron faces (codim 1) (3*4 tests)
     *       Remove duplicate points from the list
     *       Compute the convex hull polygon
     *       Return a triangulation of that polygon as intersection
     * \param geo1/geo2 The geometries to intersect
     * \param intersection Container to store the corner points of the polygon (as convex hull)
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

        const auto refElement1 = referenceElement(geo1);
        const auto refElement2 = referenceElement(geo2);

        // intersection policy for point-like intersections (used later)
        using PointPolicy = IntersectionPolicy::PointPolicy<ctype, dimworld>;

        // add intersection points of all polyhedron edges (codim dim-1) with the polygon
        for (int i = 0; i < refElement1.size(dim1-1); ++i)
        {
            const auto localEdgeGeom = refElement1.template geometry<dim1-1>(i);
            const auto p = geo1.global(localEdgeGeom.corner(0));
            const auto q = geo1.global(localEdgeGeom.corner(1));
            const auto segGeo = SegGeometry(Dune::GeometryTypes::line, std::vector<Point>{p, q});

            using PolySegTest = GeometryIntersection<Geometry2, SegGeometry, PointPolicy>;
            typename PolySegTest::Intersection polySegIntersection;
            if (PolySegTest::intersection(geo2, segGeo, polySegIntersection))
                points.emplace_back(polySegIntersection);
        }

        // add intersection points of all polygon faces (codim 1) with the polyhedron faces
        for (int i = 0; i < refElement1.size(1); ++i)
        {
            const auto faceGeo = [&]()
            {
                const auto localFaceGeo = refElement1.template geometry<1>(i);
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

            for (int j = 0; j < refElement2.size(1); ++j)
            {
                const auto localEdgeGeom = refElement2.template geometry<1>(j);
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
        intersection = grahamConvexHull<2>(points);
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
     * \param intersection Container to store the intersection result
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
 * \brief A class for polygon--polyhedron intersection in 3d space
 */
template <class Geometry1, class Geometry2, class Policy>
class GeometryIntersection<Geometry1, Geometry2, Policy, 3, 2, 3>
: public GeometryIntersection<Geometry2, Geometry1, Policy, 3, 3, 2>
{
    using Base = GeometryIntersection<Geometry2, Geometry1, Policy, 3, 3, 2>;
public:
    /*!
     * \brief Colliding polygon and convex polyhedron
     * \param geo1/geo2 The geometries to intersect
     * \param intersection If the geometries collide intersection holds the
     *        corner points of the intersection object in global coordinates.
     * \note This forwards to the polyhedron-polygon specialization with swapped arguments.
     */
    template<class P = Policy>
    static bool intersection(const Geometry1& geo1, const Geometry2& geo2, typename Base::Intersection& intersection)
    {
        return Base::intersection(geo2, geo1, intersection);
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

private:
    static constexpr ctype eps_ = 1.5e-7; // base epsilon for floating point comparisons

public:
    /*!
     *  \brief Colliding segment and convex polyhedron
     *  \note Algorithm based on the one from "Real-Time Collision Detection" by Christer Ericson,
     *        published by Morgan Kaufmann Publishers, (c) 2005 Elsevier Inc. (Chapter 5.3.6)
     * \param geo1/geo2 The geometries to intersect
     * \param intersection If the geometries collide, is holds the corner points of
     *        the intersection object in global coordinates.
     * \note This overload is used when point-like intersections are seeked
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

private:
    /*!
     * \brief Obtain local coordinates of start/end point of the intersecting segment.
     */
    template<class P = Policy, std::enable_if_t<P::dimIntersection == 1, int> = 0>
    static bool intersect_(const Geometry1& geo1, const Geometry2& geo2, ctype& tfirst, ctype& tlast)
    {
        // lambda to obtain the facet corners on geo1
        auto getFacetCorners = [] (const Geometry1& geo1)
        {
            std::vector< std::array<int, 2> > facetCorners;
            switch (geo1.corners())
            {
                case 4: // quadrilateral
                    facetCorners = {{0, 2}, {3, 1}, {1, 0}, {2, 3}};
                    break;
                case 3: // triangle
                    facetCorners = {{1, 0}, {0, 2}, {2, 1}};
                    break;
                default:
                    DUNE_THROW(Dune::NotImplemented, "Collision of segment and geometry of type "
                                   << geo1.type() << " with "<< geo1.corners() << " corners.");
            }

            return facetCorners;
        };

        const auto center1 = geo1.center();
        const auto normal1 = crossProduct(geo1.corner(1) - geo1.corner(0),
                                          geo1.corner(2) - geo1.corner(0));

        // lambda to obtain the normal on a facet
        auto computeNormal = [&center1, &normal1, &geo1] (const std::array<int, 2>& facetCorners)
        {
            const auto c0 = geo1.corner(facetCorners[0]);
            const auto c1 = geo1.corner(facetCorners[1]);
            const auto edge = c1 - c0;

            Dune::FieldVector<ctype, dimworld> n = crossProduct(edge, normal1);
            n /= n.two_norm();

            // orientation of the element is not known
            // make sure normal points outwards of element
            if ( n*(center1-c0) > 0.0 )
                n *= -1.0;

            return n;
        };

        return Detail::computeSegmentIntersection(geo1, geo2, eps_, tfirst, tlast, getFacetCorners, computeNormal);
    }

    /*!
     * \brief triangle--segment point-like intersection with points as input.
     */
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

/*!
 * \ingroup Geometry
 * \brief A class for segment--polygon intersection in 3d space
 */
template <class Geometry1, class Geometry2, class Policy>
class GeometryIntersection<Geometry1, Geometry2, Policy, 3, 1, 2>
: public GeometryIntersection<Geometry2, Geometry1, Policy, 3, 2, 1>
{
    using Base = GeometryIntersection<Geometry2, Geometry1, Policy, 3, 2, 1>;
public:
    /*!
     * \brief Colliding segment and convex polygon
     * \param geo1/geo2 The geometries to intersect
     * \param intersection If the geometries collide intersection holds the
     *        corner points of the intersection object in global coordinates.
     * \note This forwards to the polyhedron-polygon specialization with swapped arguments.
     */
    template<class P = Policy>
    static bool intersection(const Geometry1& geo1, const Geometry2& geo2, typename Base::Intersection& intersection)
    {
        return Base::intersection(geo2, geo1, intersection);
    }
};

/*!
 * \ingroup Geometry
 * \brief A class for segment--segment intersection in 3d space
 */
template <class Geometry1, class Geometry2, class Policy>
class GeometryIntersection<Geometry1, Geometry2, Policy, 3, 1, 1>
{
    enum { dimworld = 3 };
    enum { dim1 = 1 };
    enum { dim2 = 1 };

public:
    using ctype = typename Policy::ctype;
    using Point = typename Policy::Point;
    using Intersection = typename Policy::Intersection;

private:
    static constexpr ctype eps_ = 1.5e-7; // base epsilon for floating point comparisons

public:
    /*!
     * \brief Colliding two segments
     * \param geo1/geo2 The geometries to intersect
     * \param intersection Container to store corners of intersection segment
     * \note this overload is used when point-like intersections are seeked
     * \todo implement this query
     */
    template<class P = Policy, std::enable_if_t<P::dimIntersection == 0, int> = 0>
    static bool intersection(const Geometry1& geo1, const Geometry2& geo2, Intersection& intersection)
    {
        static_assert(int(dimworld) == int(Geometry2::coorddimension), "Can only collide geometries of same coordinate dimension");
        DUNE_THROW(Dune::NotImplemented, "segment-segment intersection detection for point-like intersections");
    }

    /*!
     *  \brief Colliding two segments in 3D
     * \param geo1/geo2 The geometries to intersect
     * \param intersection If the geometries collide, is holds the corner points of
     *        the intersection object in global coordinates.
     * \note This overload is used when segment-like intersections are seeked
     */
    template<class P = Policy, std::enable_if_t<P::dimIntersection == 1, int> = 0>
    static bool intersection(const Geometry1& geo1, const Geometry2& geo2, Intersection& intersection)
    {
        static_assert(int(dimworld) == int(Geometry2::coorddimension), "Can only collide geometries of same coordinate dimension");

        const auto& a = geo1.corner(0);
        const auto& b = geo1.corner(1);
        const auto ab = b-a;

        const auto& p = geo2.corner(0);
        const auto& q = geo2.corner(1);
        const auto pq = q-p;

        const auto abNorm2 = ab.two_norm2();
        const auto pqNorm2 = pq.two_norm2();

        using std::max;
        const auto eps2 = eps_*max(abNorm2, pqNorm2);

        // if the segment are not parallel there is no segment intersection
        using std::abs;
        if (crossProduct(ab, pq).two_norm2() > eps2*eps2)
            return false;

        const auto ap = (p-a);
        const auto aq = (q-a);

        // points have to be colinear
        if (crossProduct(ap, aq).two_norm2() > eps2*eps2)
            return false;

        // scaled local coordinates
        // we save the division for the normalization for now
        // and do it in the very end if we find an intersection
        auto tp = ap*ab;
        auto tq = aq*ab;

        // make sure they are sorted
        using std::swap;
        if (tp > tq)
            swap(tp, tq);

        using std::clamp;
        tp = clamp(tp, 0.0, abNorm2);
        tq = clamp(tq, 0.0, abNorm2);

        if (abs(tp-tq) < eps2)
            return false;

        intersection = Intersection({geo1.global(tp/abNorm2), geo1.global(tq/abNorm2)});
        return true;
    }
};

} // end namespace Dumux

#endif
