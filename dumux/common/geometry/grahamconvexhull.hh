/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
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
 * \ingroup Common
 * \brief A function to compute the convex hull of a point cloud
 *        and a function to triangulate the polygon area spanned by the convex hull
 */
#ifndef DUMUX_GRAHAM_CONVEX_HULL_HH
#define DUMUX_GRAHAM_CONVEX_HULL_HH

#include <vector>
#include <array>
#include <algorithm>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/math.hh>

namespace Dumux {

/*!
 * \brief Returns the orientation of a sequence a-->b-->c in one plane (defined by normal vector)
 * \return -1   if a-->b-->c forms a counter-clockwise turn (given the normal vector)
 *         +1   for a clockwise turn,
 *          0   if they are on one line (colinear)
 */
int getOrientation(const Dune::FieldVector<double, 3>& a,
                   const Dune::FieldVector<double, 3>& b,
                   const Dune::FieldVector<double, 3>& c,
                   const Dune::FieldVector<double, 3>& normal)
{
    const auto d = b-a;
    const auto e = c-b;
    const auto f = Dumux::crossProduct(d, e);
    const auto area = f*normal;
    return Dumux::sign(-area);
}

/*!
 * \brief Compute the points making up the convex hull around the given set of unordered points
 * \note We assume that all points are coplanar
 * \note this algorithm changes the order of the given points a bit
 *       as they are unordered anyway this shouldn't matter too much
 */
std::vector<Dune::FieldVector<double, 3>>
grahamConvexHull2d3d(std::vector<Dune::FieldVector<double, 3>>& points)
{
    using Point = Dune::FieldVector<double, 3>;
    std::vector<Point> convexHull; convexHull.reserve(50);

    // return empty convex hull
    if (points.size() < 3)
        return convexHull;

    // return the points (already just one triangle)
    if (points.size() == 3)
        return points;

    // get the normal vector of the plane
    const auto a = points[1] - points[0];
    const auto b = points[2] - points[0];
    auto normal = Dumux::crossProduct(a, b);
    normal /= normal.two_norm();

    // find the element with the smallest y coordinate (if y is the same, smallest x coordinate)
    auto minIt = std::min_element(points.begin(), points.end(),
                                  [](const auto& a, const auto& b)
                                  { return a[1] != b[1] ? a[1] < b[1] : a[0] < b[0]; });

    // swap the smallest element to the front
    std::iter_swap(minIt, points.begin());

    // choose the first (min element) as the pivot point
    // sort in counter-clockwise order around pivot point
    const auto pivot = points[0];
    std::sort(points.begin()+1, points.end(), [&](const auto& a, const auto& b)
              {
                  int order = getOrientation(pivot, a, b, normal);
                  if (order == 0)
                       return (a-pivot).two_norm() < (b-pivot).two_norm();
                  else
                       return (order == -1);
              });

    // push the first three points
    convexHull.push_back(points[0]);
    convexHull.push_back(points[1]);
    convexHull.push_back(points[2]);

    // This is the heart of the algorithm
    // Pop_back until the last point in the queue forms a counter-clockwise oriented line
    // with the new vertices. Then add new points to the queue.
    for (std::size_t i = 3; i < points.size(); ++i)
    {
        Point p = convexHull.back();
        convexHull.pop_back();
        // keep popping until the orientation a->b->currentp is counter-clockwise
        while (getOrientation(convexHull.back(), p, points[i], normal) != -1)
        {
            p = convexHull.back();
            convexHull.pop_back();
        }

        // add back the last popped point and this point
        convexHull.emplace_back(std::move(p));
        convexHull.push_back(points[i]);
    }

    return convexHull;
}

/*!
 * \brief Triangulate area given points of the convex hull
 * \note Assumes all points of the convex hull are coplanar
 * \note This inserts a mid point and connects all corners with that point to triangles
 */
std::vector<std::array<Dune::FieldVector<double, 3>, 3> >
triangulateConvexHull(const std::vector<Dune::FieldVector<double, 3>>& convexHull)
{
    using Point = Dune::FieldVector<double, 3>;
    using Triangle = std::array<Point, 3>;

    if (convexHull.size() < 3)
        DUNE_THROW(Dune::InvalidStateException, "Try to triangulate point cloud with less than 3 points!");

    if (convexHull.size() == 3)
        return std::vector<Triangle>(1, {convexHull[0], convexHull[1], convexHull[2]});

    Point midPoint(0.0);
    for (const auto p : convexHull)
        midPoint += p;
    midPoint /= convexHull.size();

    std::vector<Triangle> triangulation;
    triangulation.reserve(convexHull.size());

    for (std::size_t i = 0; i < convexHull.size()-1; ++i)
        triangulation.emplace_back(Triangle{midPoint, convexHull[i], convexHull[i+1]});

    triangulation.emplace_back(Triangle{midPoint, convexHull[convexHull.size()-1], convexHull[0]});

    return triangulation;
}

} // end namespace Dumux

# endif
