// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Geometry
 * \brief Functionality to triangulate point clouds
 * \note Most of the implemented algorithms currently do not scale for large point clouds
 *       and are only meant to be used with a small number of points
 */
#ifndef DUMUX_GEOMETRY_TRIANGULATION_HH
#define DUMUX_GEOMETRY_TRIANGULATION_HH

#include <vector>
#include <array>
#include <algorithm>
#include <iterator>
#include <type_traits>
#include <tuple>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/math.hh>
#include <dumux/geometry/intersectspointsimplex.hh>
#include <dumux/geometry/grahamconvexhull.hh>

namespace Dumux {
namespace TriangulationPolicy {

//! Policy that expects a point cloud that represents a convex
//! hull. Inserts the mid point and connects it to the points of the hull.
struct MidPointPolicy {};

//! Policy that first finds the convex hull
//! and then uses the mid point policy to construct a triangulation
struct ConvexHullPolicy {};

//! Delaunay-type triangulations.
struct DelaunayPolicy {};

#ifndef DOXYGEN
namespace Detail {
using DefaultDimPolicies = std::tuple<DelaunayPolicy, MidPointPolicy, ConvexHullPolicy>;
} // end namespace Detail
#endif

//! Default policy for a given dimension
template<int dim, int dimWorld>
using DefaultPolicy = std::tuple_element_t<dim-1, Detail::DefaultDimPolicies>;

} // end namespace TriangulationPolicy

/*!
 * \ingroup Geometry
 * \brief The default data type to store triangulations
 * \note Stores each simplex separate and without connectivity information
 * \note We usually use this for sub-triangulation to allow for quadrature rules and interpolation on intersection geometries
 *       This is neither meant to be used to store large amounts of data nor as a mesh-like object
 */
template<int dim, int dimWorld, class ctype>
using Triangulation = std::vector< std::array<Dune::FieldVector<ctype, dimWorld>, dim+1> >;

/*!
 * \ingroup Geometry
 * \brief Triangulate area given points of a convex hull (1d)
 * \tparam dim Specifies the dimension of the resulting triangulation
 * \tparam dimWorld The dimension of the coordinate space
 * \tparam Policy Specifies the algorithm to be used for triangulation
 *
 * \note this specialization is for 1d discretization using segments
 * \note sorts the points and connects them via segments
 */
template< int dim, int dimWorld, class Policy = TriangulationPolicy::DefaultPolicy<dim, dimWorld>,
          class RandomAccessContainer,
          std::enable_if_t< std::is_same_v<Policy, TriangulationPolicy::DelaunayPolicy>
                            && dim == 1, int> = 0 >
inline Triangulation<dim, dimWorld, typename RandomAccessContainer::value_type::value_type>
triangulate(const RandomAccessContainer& points)
{
    using ctype = typename RandomAccessContainer::value_type::value_type;
    using Point = Dune::FieldVector<ctype, dimWorld>;

    static_assert(std::is_same_v<typename RandomAccessContainer::value_type, Point>,
                  "Triangulation expects Dune::FieldVector as point type");

    if (points.size() == 2)
        return Triangulation<dim, dimWorld, ctype>({ {points[0], points[1]} });

    //! \todo sort points and create polyline
    assert(points.size() > 1);
    DUNE_THROW(Dune::NotImplemented, "1d triangulation for point cloud size > 2");
}

/*!
 * \ingroup Geometry
 * \brief Triangulate area given points of a convex hull (2d)
 * \tparam dim Specifies the dimension of the resulting triangulation
 * \tparam dimWorld The dimension of the coordinate space
 * \tparam Policy Specifies the algorithm to be used for triangulation
 *
 * \note this specialization is for 2d triangulations using mid point policy
 * \note Assumes all points of the convex hull are coplanar
 * \note This inserts a mid point and connects all corners with that point to triangles
 * \note Assumes points are given as a ordered sequence representing the polyline forming the convex hull
 */
template< int dim, int dimWorld, class Policy = TriangulationPolicy::DefaultPolicy<dim, dimWorld>,
          class RandomAccessContainer,
          std::enable_if_t< std::is_same_v<Policy, TriangulationPolicy::MidPointPolicy>
                            && dim == 2, int> = 0 >
inline Triangulation<dim, dimWorld, typename RandomAccessContainer::value_type::value_type>
triangulate(const RandomAccessContainer& convexHullPoints)
{
    using ctype = typename RandomAccessContainer::value_type::value_type;
    using Point = Dune::FieldVector<ctype, dimWorld>;
    using Triangle = std::array<Point, 3>;

    static_assert(std::is_same_v<typename RandomAccessContainer::value_type, Point>,
                  "Triangulation expects Dune::FieldVector as point type");

    if (convexHullPoints.size() < 3)
        DUNE_THROW(Dune::InvalidStateException, "Try to triangulate point cloud with less than 3 points!");

    if (convexHullPoints.size() == 3)
        return std::vector<Triangle>(1, {convexHullPoints[0], convexHullPoints[1], convexHullPoints[2]});

    Point midPoint(0.0);
    for (const auto& p : convexHullPoints)
        midPoint += p;
    midPoint /= convexHullPoints.size();

    std::vector<Triangle> triangulation;
    triangulation.reserve(convexHullPoints.size());

    for (std::size_t i = 0; i < convexHullPoints.size()-1; ++i)
        triangulation.emplace_back(Triangle{midPoint, convexHullPoints[i], convexHullPoints[i+1]});

    triangulation.emplace_back(Triangle{midPoint, convexHullPoints[convexHullPoints.size()-1], convexHullPoints[0]});

    return triangulation;
}

/*!
 * \ingroup Geometry
 * \brief Triangulate area given points (2d)
 * \tparam dim Specifies the dimension of the resulting triangulation
 * \tparam dimWorld The dimension of the coordinate space
 * \tparam Policy Specifies the algorithm to be used for triangulation
 *
 * \note this specialization is for 2d triangulations using the convex hull policy
 *       That means we first construct the convex hull of the points and then
 *       triangulate the convex hull using the midpoint policy
 * \note Assumes points are unique and not all colinear (will throw a Dune::InvalidStateException)
 */
template< int dim, int dimWorld, class Policy = TriangulationPolicy::DefaultPolicy<dim, dimWorld>,
          class RandomAccessContainer,
          std::enable_if_t< std::is_same_v<Policy, TriangulationPolicy::ConvexHullPolicy>
                            && dim == 2, int> = 0 >
inline Triangulation<dim, dimWorld, typename RandomAccessContainer::value_type::value_type>
triangulate(const RandomAccessContainer& points)
{
    const auto convexHullPoints = grahamConvexHull<2>(points);
    return triangulate<dim, dimWorld, TriangulationPolicy::MidPointPolicy>(convexHullPoints);
}

/*!
 * \ingroup Geometry
 * \brief Triangulate volume given a point cloud (3d)
 * \tparam dim Specifies the dimension of the resulting triangulation
 * \tparam dimWorld The dimension of the coordinate space
 * \tparam Policy Specifies the algorithm to be used for triangulation
 *
 * \note this specialization is for 3d triangulations using the convex hull policy
 * \note Assumes points are unique and not all coplanar
 */
template< int dim, int dimWorld, class Policy = TriangulationPolicy::DefaultPolicy<dim, dimWorld>,
          class RandomAccessContainer,
          std::enable_if_t< std::is_same_v<Policy, TriangulationPolicy::ConvexHullPolicy>
                            && dim == 3, int> = 0 >
inline Triangulation<dim, dimWorld, typename RandomAccessContainer::value_type::value_type>
triangulate(const RandomAccessContainer& points)
{
    using ctype = typename RandomAccessContainer::value_type::value_type;
    using Point = Dune::FieldVector<ctype, dimWorld>;
    using Tetrahedron = std::array<Point, 4>;

    static_assert(std::is_same_v<typename RandomAccessContainer::value_type, Point>,
                  "Triangulation expects Dune::FieldVector as point type");

    const auto numPoints = points.size();
    if (numPoints < 4)
        DUNE_THROW(Dune::InvalidStateException, "Trying to create 3D triangulation of point cloud with less than 4 points!");

    if (numPoints == 4)
        return std::vector<Tetrahedron>(1, {points[0], points[1], points[2], points[3]});

    // compute the mid point of the point cloud (not the midpoint of the convex hull but this
    // should not matter too much for the applications we have in mind here)
    Point midPoint(0.0);
    Point lowerLeft(1e100);
    Point upperRight(-1e100);
    for (const auto& p : points)
    {
        midPoint += p;
        for (int i = 0; i < dimWorld; ++i)
        {
            using std::max; using std::min;
            lowerLeft[i] = min(p[i], lowerLeft[i]);
            upperRight[i] = max(p[i], upperRight[i]);
        }
    }
    midPoint /= numPoints;

    auto magnitude = 0.0;
    using std::max;
    for (int i = 0; i < dimWorld; ++i)
        magnitude = max(upperRight[i] - lowerLeft[i], magnitude);
    const auto eps = 1e-7*magnitude;
    const auto eps2 = eps*eps;
    const auto epsDist = 1e-8*magnitude;

    // faces of the convex hull given by the sorted indices of their (coplanar) points
    std::vector<std::vector<int>> faces;
    faces.reserve(2*numPoints);
    std::vector<int> face;
    face.reserve(numPoints);

    const auto isOnFace = [&](int i, int j, int k)
    {
        return std::any_of(faces.begin(), faces.end(), [&](const auto& f)
        {
            return std::binary_search(f.begin(), f.end(), i)
                && std::binary_search(f.begin(), f.end(), j)
                && std::binary_search(f.begin(), f.end(), k);
        });
    };

    // Two faces of a convex polytope share at most an edge, so two faces sharing
    // three points that are not on a line are the same face, found from triangles
    // that classified coplanar points differently. Points on a common edge deviate
    // from it by about the coplanarity tolerance epsDist, far below the threshold used here.
    const auto isSameFace = [&](const auto& f, const auto& g)
    {
        std::vector<int> shared;
        std::set_intersection(f.begin(), f.end(), g.begin(), g.end(), std::back_inserter(shared));
        if (shared.size() < 3)
            return false;

        std::size_t a = 0, b = 1;
        for (std::size_t m = 0; m < shared.size(); ++m)
            for (std::size_t n = m+1; n < shared.size(); ++n)
                if ((points[shared[n]] - points[shared[m]]).two_norm2() > (points[shared[b]] - points[shared[a]]).two_norm2())
                { a = m; b = n; }

        const auto line = points[shared[b]] - points[shared[a]];
        const auto lineLength = line.two_norm();
        return std::any_of(shared.begin(), shared.end(), [&](int c)
        { return crossProduct(line, points[c] - points[shared[a]]).two_norm() > 1e-6*magnitude*lineLength; });
    };

    // brute force algorithm: Try all possible triangles and check
    // if they are triangles of the convex hull. This is achieved by checking if all
    // other points are in the same half-space (side)
    // the algorithm is O(n^4), so only do this for very small point clouds
    for (int i = 0; i < numPoints; ++i)
    {
        for (int j = i+1; j < numPoints; ++j)
        {
            for (int k = j+1; k < numPoints; ++k)
            {
                // this triangle lies on a face that has already been found
                if (isOnFace(i, j, k))
                    continue;

                const auto pointI = points[i];
                const auto ab = points[j] - pointI;
                const auto ac = points[k] - pointI;
                const auto normal = crossProduct(ab, ac);

                // check if points are colinear and we can't form a triangle
                // if so, skip this triangle
                if (normal.two_norm2() < eps2*eps2)
                    continue;

                // check if this triangle ijk is admissible which means
                // it is on the convex hull (all other points in the cloud are in the same half-space/side)
                // and collect the points that are coplanar with it
                face.assign({ i, j, k });
                const bool isAdmissible = [&]()
                {
                    int marker = 0; // 0 means undecided side (or coplanar)
                    for (int m = 0; m < numPoints; ++m)
                    {
                        if (m != i && m != j && m != k)
                        {
                            // check scalar product with surface normal to decide side
                            const auto ad = points[m] - pointI;
                            const auto sp = normal*ad;

                            // if the sign changes wrt the previous sign, the triangle is not part of the convex hull
                            // sp is a triple product, so divide by the normal to obtain a distance
                            using std::abs; using std::signbit;
                            const bool coplanar = abs(sp) < epsDist*normal.two_norm();
                            int newMarker = coplanar ? 0 : signbit(sp) ? -1 : 1;

                            // make decision for a side as soon as the next marker is != 0
                            // keep track of the previous marker
                            if (marker == 0 && newMarker != 0)
                                marker = newMarker;

                            // if marker flips, not all points are on one side
                            // zero marker (undecided side / coplanar point) shouldn't abort the process
                            if (newMarker != 0 && marker != newMarker)
                                return false;

                            if (coplanar)
                                face.push_back(m);
                        }
                    }

                    // we require that not all points are coplanar, so
                    // there will be always at least one non-coplanar point
                    // to check on which side the other points are.
                    // Hence, once we get here, the triangle or coplanar point cluster is part of the convex hull
                    return true;
                }();

                if (isAdmissible)
                {
                    std::sort(face.begin(), face.end());
                    faces.push_back(face);
                }
            }
        }
    }

    for (bool merged = true; merged; )
    {
        merged = false;
        for (std::size_t f = 0; f < faces.size() && !merged; ++f)
        {
            for (std::size_t g = f+1; g < faces.size() && !merged; ++g)
            {
                if (isSameFace(faces[f], faces[g]))
                {
                    std::vector<int> united;
                    std::set_union(faces[f].begin(), faces[f].end(), faces[g].begin(), faces[g].end(), std::back_inserter(united));
                    faces[f] = std::move(united);
                    faces.erase(faces.begin() + g);
                    merged = true;
                }
            }
        }
    }

    // reserve memory conservatively to avoid reallocation
    std::vector<Tetrahedron> triangulation;
    triangulation.reserve(numPoints);

    // connect the triangles of each face to the mid point
    const auto addTetrahedron = [&](const Point& a, const Point& b, const Point& c)
    {
        const auto normal = crossProduct(b - a, c - a);
        using std::signbit;
        if (signbit(normal*(midPoint - a)))
            triangulation.emplace_back(Tetrahedron{ a, c, b, midPoint });
        else
            triangulation.emplace_back(Tetrahedron{ a, b, c, midPoint });
    };

    std::vector<Point> facePoints;
    for (const auto& f : faces)
    {
        if (f.size() == 3)
            addTetrahedron(points[f[0]], points[f[1]], points[f[2]]);

        // a cluster of coplanar points forms one of the faces of the convex hull:
        // compute its (2d) convex hull first and triangulate
        else
        {
            facePoints.clear();
            for (int idx : f)
                facePoints.push_back(points[idx]);
            for (const auto& triangle : triangulate<2, 3, TriangulationPolicy::ConvexHullPolicy>(facePoints))
                addTetrahedron(triangle[0], triangle[1], triangle[2]);
        }
    }

    // sanity check: if points are not coplanar, then using the mid point policy, we get at least 4 tetrahedrons
    if (triangulation.size() < 4)
        DUNE_THROW(Dune::InvalidStateException, "Something went wrong with the triangulation!");

    return triangulation;
}

} // end namespace Dumux

#endif
