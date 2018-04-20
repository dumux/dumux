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
 * \ingroup Common
 * \brief Create Dune geometries from user-specified points
 */
#ifndef DUMUX_MAKE_GEOMETRY_HH
#define DUMUX_MAKE_GEOMETRY_HH

#include <vector>
#include <array>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/version.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dumux/common/math.hh>
#include <dumux/common/geometry/grahamconvexhull.hh>

namespace Dumux {

//! Checks if four points lie within the same plane.
template<class Scalar>
bool pointsAreCoplanar(const std::vector<Dune::FieldVector<Scalar, 3>>& points, Scalar eps = 1e-20)
{
    assert(points.size() == 4);
    // (see "Real-Time Collision Detection" by Christer Ericson)
    Dune::FieldMatrix<Scalar, 4, 4> M;
    for(int i = 0; i < 3; ++i )
        M[i] = {points[0][i], points[1][i], points[2][i], points[3][i]};
    M[3] = {1.0, 1.0, 1.0, 1.0};

    using std::abs;
    return abs(M.determinant()) < eps;
}

/*!
 * \brief Returns a vector of points following the dune ordering.
 *        Convenience method that creates a temporary object in case no array of orientations is desired.
 *
 * \param points The user-specified vector of points (potentially in wrong order).
 */
template<class Scalar>
std::vector<Dune::FieldVector<Scalar, 3>> getReorderedPoints(const std::vector<Dune::FieldVector<Scalar, 3>>& points)
{
    std::array<int, 4> tmp;
    return getReorderedPoints(points, tmp);
}

/*!
 * \brief Returns a vector of points following the dune ordering.
 *
 * \param points The user-specified vector of points (potentially in wrong order).
 * \param orientations An array of orientations that can be useful for further processing.
 */
template<class Scalar>
std::vector<Dune::FieldVector<Scalar, 3>> getReorderedPoints(const std::vector<Dune::FieldVector<Scalar, 3>>& points,
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
        using GlobalPosition = Dune::FieldVector<Scalar, 3>;
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
 * \brief Creates a dune quadrilateral geometry given 4 corner points.
 *
 * \tparam Scalar The Scalar type.
 * \tparam enableSanityCheck Turn on/off sanity check and reordering of points
 * \param points The user-specified vector of points (potentially in wrong order).
 */
template<class Scalar, bool enableSanityCheck = true>
auto makeDuneQuadrilaterial(const std::vector<Dune::FieldVector<Scalar, 3>>& points)
{
    assert(points.size() == 4 && "A quadrilateral needs 4 corner points!");

    using GlobalPosition = Dune::FieldVector<Scalar, 3>;
    static constexpr auto coordDim = GlobalPosition::dimension;
    static constexpr auto dim = coordDim-1;
    using GeometryType = Dune::MultiLinearGeometry<Scalar, dim, coordDim>;

    // if no sanity check if desired, use the given points directly to construct the geometry
    if(!enableSanityCheck)
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,6)
        return GeometryType(Dune::GeometryTypes::quadrilateral, points);
#else
    {
        static Dune::GeometryType gt(Dune::GeometryType::cube, dim);
        return GeometryType(gt, points);
    }
#endif

    // otherwise, perform a number of checks and corrections
    if(!pointsAreCoplanar(points))
        DUNE_THROW(Dune::InvalidStateException, "Points do not lie within a plane");

    // make sure points conform with dune ordering
    std::array<int, 4> orientations;
    auto corners = getReorderedPoints(points, orientations);

    if(std::any_of(orientations.begin(), orientations.end(), [](auto i){ return i == 0; }))
        DUNE_THROW(Dune::InvalidStateException, "More than two points lie on the same line.");

#if DUNE_VERSION_NEWER(DUNE_COMMON,2,6)
    const auto quadrilateral = GeometryType(Dune::GeometryTypes::quadrilateral, corners);
#else
    static Dune::GeometryType gt(Dune::GeometryType::cube, dim);
    const auto quadrilateral = GeometryType(gt, corners);
#endif

    const auto eps = 1e-20;
    if(quadrilateral.volume() < eps)
        DUNE_THROW(Dune::InvalidStateException, "Something went wrong, geometry has area of zero");

    return quadrilateral;
}



} // end namespace Dumux

#endif // DUMUX_MAKE_GEOMETRY_HH
