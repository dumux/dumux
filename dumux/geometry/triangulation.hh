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
 * \brief Functionality to triangulate point clouds
 */
#ifndef DUMUX_GEOMETRY_TRIANGULATION_HH
#define DUMUX_GEOMETRY_TRIANGULATION_HH

#include <vector>
#include <array>
#include <algorithm>
#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/math.hh>

namespace Dumux {
namespace TriangulationPolicy {

//! Policy that expects a point cloud that represents a convex
//! hull. Inserts the mid point and connects it to the points of the hull.
struct MidPointPolicy {};

//! Delaunay-type triangulations.
struct DelaunayPolicy {};

template<int dim, int dimWorld>
using DefaultPolicy = std::conditional_t< dim >= 2, MidPointPolicy, DelaunayPolicy >;

} // end namespace TriangulationPolicy

//! The data type to store triangulations
template<int dim, int dimWorld, class ctype>
using Triangulation = std::vector< std::array<Dune::FieldVector<ctype, dimWorld>, dim+1> >;

/*!
 * \ingroup Geometry
 * \brief Triangulate area given points of a convex hull
 * \tparam dim Specifies the dimension of the resulting triangulation (2 or 3)
 * \tparam dimWorld The dimension of the coordinate space
 * \tparam Policy Specifies the algorithm to be used for triangulation
 *
 * \note this specialization is for 2d triangulations using mid point policy
 * \note Assumes all points of the convex hull are coplanar
 * \note This inserts a mid point and connects all corners with that point to triangles
 */
template< int dim, int dimWorld, class Policy = TriangulationPolicy::DefaultPolicy<dim, dimWorld>,
          class RandomAccessContainer,
          std::enable_if_t< std::is_same<Policy, TriangulationPolicy::MidPointPolicy>::value
                            && dim == 2, int> = 0 >
inline Triangulation<dim, dimWorld, typename RandomAccessContainer::value_type::value_type>
triangulate(const RandomAccessContainer& convexHull)
{
    using ctype = typename RandomAccessContainer::value_type::value_type;
    using Point = Dune::FieldVector<ctype, dimWorld>;
    using Triangle = std::array<Point, 3>;

    static_assert(std::is_same<typename RandomAccessContainer::value_type, Point>::value,
                  "Triangulation expects Dune::FieldVector as point type");

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

/*!
 * \ingroup Geometry
 * \brief Triangulate area given points of a convex hull
 * \tparam dim Specifies the dimension of the resulting triangulation (2 or 3)
 * \tparam dimWorld The dimension of the coordinate space
 * \tparam Policy Specifies the algorithm to be used for triangulation
 *
 * \note this specialization is for 1d discretization using segments
 * \note sorts the points and connects them via segments
 */
template< int dim, int dimWorld, class Policy = TriangulationPolicy::DefaultPolicy<dim, dimWorld>,
          class RandomAccessContainer,
          std::enable_if_t< std::is_same<Policy, TriangulationPolicy::DelaunayPolicy>::value
                            && dim == 1, int> = 0 >
inline Triangulation<dim, dimWorld, typename RandomAccessContainer::value_type::value_type>
triangulate(const RandomAccessContainer& points)
{
    using ctype = typename RandomAccessContainer::value_type::value_type;
    using Point = Dune::FieldVector<ctype, dimWorld>;

    static_assert(std::is_same<typename RandomAccessContainer::value_type, Point>::value,
                  "Triangulation expects Dune::FieldVector as point type");

    if (points.size() == 2)
        return Triangulation<dim, dimWorld, ctype>({ {points[0], points[1]} });

    //! \todo sort points and create polyline
    assert(points.size() > 1);
    DUNE_THROW(Dune::NotImplemented, "1d triangulation for point cloud size > 2");
}

} // end namespace Dumux

#endif
