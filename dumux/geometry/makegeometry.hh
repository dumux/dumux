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
 * \brief Create Dune geometries from user-specified points
 */
#ifndef DUMUX_GEOMETRY_MAKE_GEOMETRY_HH
#define DUMUX_GEOMETRY_MAKE_GEOMETRY_HH

#include <vector>
#include <array>
#include <limits>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/exceptions.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dumux/common/math.hh>
#include <dumux/geometry/grahamconvexhull.hh>

namespace Dumux {

/*!
 * \ingroup Geometry
 * \brief Checks if four points lie within the same plane
 */
template<class CoordScalar>
bool pointsAreCoplanar(const std::vector<Dune::FieldVector<CoordScalar, 3>>& points, const CoordScalar scale)
{
    if (points.size() != 4)
        DUNE_THROW(Dune::InvalidStateException, "Check only works for 4 points!");

    // (see "Real-Time Collision Detection" by Christer Ericson)
    Dune::FieldMatrix<CoordScalar, 4, 4> M;
    for (int i = 0; i < 3; ++i )
        M[i] = {points[0][i], points[1][i], points[2][i], points[3][i]};
    M[3] = {1.0*scale, 1.0*scale, 1.0*scale, 1.0*scale};

    using std::abs;
    return abs(M.determinant()) < 1.5e-7*scale*scale*scale*scale;
}

/*!
 * \ingroup Geometry
 * \brief  Checks if four points lie within the same plane.
 */
template<class CoordScalar>
bool pointsAreCoplanar(const std::vector<Dune::FieldVector<CoordScalar, 3>>& points)
{
    Dune::FieldVector<CoordScalar, 3> bBoxMin(std::numeric_limits<CoordScalar>::max());
    Dune::FieldVector<CoordScalar, 3> bBoxMax(std::numeric_limits<CoordScalar>::lowest());
    for (const auto& p : points)
    {
        for (int i=0; i<3; i++)
        {
            using std::min;
            using std::max;
            bBoxMin[i] = min(bBoxMin[i], p[i]);
            bBoxMax[i] = max(bBoxMax[i], p[i]);
        }
    }

    const auto size = (bBoxMax - bBoxMin).two_norm();

    return pointsAreCoplanar(points, size);
}

/*!
 * \ingroup Geometry
 * \brief Returns a vector of points following the dune ordering.
 *        Convenience method that creates a temporary object in case no array of orientations is desired.
 *
 * \param points The user-specified vector of points (potentially in wrong order).
 */
template<class CoordScalar>
std::vector<Dune::FieldVector<CoordScalar, 3>> getReorderedPoints(const std::vector<Dune::FieldVector<CoordScalar, 3>>& points)
{
    std::array<int, 4> tmp;
    return getReorderedPoints(points, tmp);
}

/*!
 * \ingroup Geometry
 * \brief Returns a vector of points following the dune ordering.
 *
 * \param points The user-specified vector of points (potentially in wrong order).
 * \param orientations An array of orientations that can be useful for further processing.
 */
template<class CoordScalar>
std::vector<Dune::FieldVector<CoordScalar, 3>> getReorderedPoints(const std::vector<Dune::FieldVector<CoordScalar, 3>>& points,
                                                             std::array<int, 4>& orientations)
{
    if(points.size() == 4)
    {
        auto& p0 = points[0];
        auto& p1 = points[1];
        auto& p2 = points[2];
        auto& p3 = points[3];

        // check if the points define a proper quadrilateral
        const auto normal = crossProduct((p1 - p0), (p2 - p0));

        orientations = { getOrientation(p0, p3, p2, normal),
                         getOrientation(p0, p3, p1, normal),
                         getOrientation(p2, p1, p0, normal),
                         getOrientation(p2, p1, p3, normal) };


        // check if the points follow the dune ordering (see http://www.dcs.gla.ac.uk/~pat/52233/slides/Geometry1x1.pdf)
        const bool diagonalsIntersect = (orientations[0] != orientations[1]) && (orientations[2] != orientations[3]);

        // the points conform with the dune ordering
        if(diagonalsIntersect)
            return points;

        // the points do not conform with the dune ordering, re-order
        using GlobalPosition = Dune::FieldVector<CoordScalar, 3>;
        if(!diagonalsIntersect && orientations[0] == 1)
            return std::vector<GlobalPosition>{p1, p0, p2, p3};
        else if(!diagonalsIntersect && orientations[0] == -1)
            return std::vector<GlobalPosition>{p3, p1, p0, p2};
        else
            DUNE_THROW(Dune::InvalidStateException, "Could not reorder points");
    }
    else
        DUNE_THROW(Dune::NotImplemented, "Reorder for " << points.size() << " points.");
}

/*!
 * \ingroup Geometry
 * \brief Creates a dune quadrilateral geometry given 4 corner points.
 *
 * \tparam CoordScalar The CoordScalar type.
 * \tparam enableSanityCheck Turn on/off sanity check and reordering of points
 * \param points The user-specified vector of points (potentially in wrong order).
 */
template<class CoordScalar, bool enableSanityCheck = true>
auto makeDuneQuadrilaterial(const std::vector<Dune::FieldVector<CoordScalar, 3>>& points)
{
    if (points.size() != 4)
        DUNE_THROW(Dune::InvalidStateException, "A quadrilateral needs 4 corner points!");

    using GlobalPosition = Dune::FieldVector<CoordScalar, 3>;
    static constexpr auto coordDim = GlobalPosition::dimension;
    static constexpr auto dim = coordDim-1;
    using GeometryType = Dune::MultiLinearGeometry<CoordScalar, dim, coordDim>;

    // if no sanity check if desired, use the given points directly to construct the geometry
    if (!enableSanityCheck)
        return GeometryType(Dune::GeometryTypes::quadrilateral, points);

    // compute size
    Dune::FieldVector<CoordScalar, 3> bBoxMin(std::numeric_limits<CoordScalar>::max());
    Dune::FieldVector<CoordScalar, 3> bBoxMax(std::numeric_limits<CoordScalar>::lowest());
    for (const auto& p : points)
    {
        for (int i = 0; i < 3; i++)
        {
            using std::min;
            using std::max;
            bBoxMin[i] = min(bBoxMin[i], p[i]);
            bBoxMax[i] = max(bBoxMax[i], p[i]);
        }
    }

    const auto size = (bBoxMax - bBoxMin).two_norm();

    // otherwise, perform a number of checks and corrections
    if (!pointsAreCoplanar(points, size))
        DUNE_THROW(Dune::InvalidStateException, "Points do not lie within a plane");

    auto corners = grahamConvexHull<2>(points);
    if (corners.size() != 4)
        DUNE_THROW(Dune::InvalidStateException, "Points do not span a strictly convex polygon!");

    // make sure points conform with dune ordering
    std::array<int, 4> orientations;
    corners = getReorderedPoints(corners, orientations);

    if (std::any_of(orientations.begin(), orientations.end(), [](auto i){ return i == 0; }))
        DUNE_THROW(Dune::InvalidStateException, "More than two points lie on the same line.");

    const auto quadrilateral = GeometryType(Dune::GeometryTypes::quadrilateral, corners);

    const auto eps = 1e-7;
    if (quadrilateral.volume() < eps*size*size)
        DUNE_THROW(Dune::InvalidStateException, "Something went wrong, geometry has area of zero");

    return quadrilateral;
}



} // end namespace Dumux

#endif
