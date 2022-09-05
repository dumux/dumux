// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
 * \brief Compute the volume of several common geometry types
 */
#ifndef DUMUX_GEOMETRY_VOLUME_HH
#define DUMUX_GEOMETRY_VOLUME_HH

#include <cmath>
#include <limits>
#include <type_traits>
#include <ranges>

#include <dune/common/exceptions.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dumux/common/math.hh>

namespace Dumux {

// TODO: properly define the concept of a vector (subset of Dune::FieldVector)
//       that we rely upon here!
// TODO: maybe provide an overload with the geometry id as template argument
//       such that one can avoid the ifs in case the type is known at compile-time

/*!
 * \ingroup Geometry
 * \brief Compute the volume of several common geometry types
 * \param type the geometry type
 * \param c a function returning the ith corner (in Dune reference element order)
 *         e.g. `[&](unsigned int i){ return corners[i]; }`, where `corners` is a
 *         random access container storing the corners, and the returned corner is stored
 *         in a container (e.g. Dune::FieldVector) that has `::value_type` and `::dimension`.
 * \tparam dim the dimension of the geometry
 * \tparam CornerF the function type (is deduced)
 * \return volume of the geometry or NaN signalling not implemented
 * \note This is only correct for convex polytopes (flat sides)
 */
template<int dim, std::invocable<unsigned int> CornerF>
auto convexPolytopeVolume(Dune::GeometryType type, const CornerF& c)
{
    using ctype = typename std::decay_t<decltype(c(0))>::value_type;
    static constexpr int coordDim = std::decay_t<decltype(c(0))>::dimension;
    static_assert(coordDim >= dim, "Coordinate dimension has to be larger than geometry dimension");

    // not implemented for coordinate dimension larger than 3
    if constexpr (coordDim > 3)
        return std::numeric_limits<ctype>::quiet_NaN();

    if constexpr (dim == 0)
        return 1.0;

    else if constexpr (dim == 1)
        return (c(1)-c(0)).two_norm();

    else if constexpr (dim == 2)
    {
        if (type == Dune::GeometryTypes::triangle)
        {
            if constexpr (coordDim == 2)
            {
                // make sure we are using positive volumes
                // the cross product of edge vectors might be negative,
                // depending on the element orientation
                using std::abs;
                return 0.5*abs(Dumux::crossProduct(c(1)-c(0), c(2)-c(0)));
            }
            else // coordDim == 3
                return 0.5*Dumux::crossProduct(c(1)-c(0), c(2)-c(0)).two_norm();

        }
        else if (type == Dune::GeometryTypes::quadrilateral)
        {
            if constexpr (coordDim == 2)
            {
                // make sure we are using positive volumes
                // the cross product of diagonals might be negative,
                // depending on the element orientation
                using std::abs;
                return 0.5*abs(Dumux::crossProduct(c(3)-c(0), c(2)-c(1)));
            }
            else // coordDim == 3
                return 0.5*Dumux::crossProduct(c(3)-c(0), c(2)-c(1)).two_norm();

        }
        else
            return std::numeric_limits<ctype>::quiet_NaN();
    }

    else if constexpr (dim == 3)
    {
        if (type == Dune::GeometryTypes::tetrahedron)
        {
            using std::abs;
            return 1.0/6.0 * abs(
                Dumux::tripleProduct(c(3)-c(0), c(1)-c(0), c(2)-c(0))
            );
        }
        else if (type == Dune::GeometryTypes::hexahedron)
        {
            // after Grandy 1997, Efficient computation of volume of hexahedron
            const auto v = c(7)-c(0);
            using std::abs;
            return 1.0/6.0 * (
                abs(Dumux::tripleProduct(v, c(1)-c(0), c(3)-c(5)))
                + abs(Dumux::tripleProduct(v, c(4)-c(0), c(5)-c(6)))
                + abs(Dumux::tripleProduct(v, c(2)-c(0), c(6)-c(3)))
            );
        }
        else if (type == Dune::GeometryTypes::pyramid)
        {
            // 1/3 * base * height
            // for base see case Dune::GeometryTypes::quadrilateral above
            // = 1/3 * (1/2 * norm(ADxBC)) * ((ADxBC)/norm(AD x BC) ⋅ AE)
            // = 1/6 * (AD x BC) ⋅ AE
            using std::abs;
            return 1.0/6.0 * abs(
                Dumux::tripleProduct(c(3)-c(0), c(2)-c(1), c(4)-c(0))
            );
        }
        else if (type == Dune::GeometryTypes::prism)
        {
            // compute as sum of a pyramid (0-1-3-4-5) and a tetrahedron (2-0-1-5)
            using std::abs;
            return 1.0/6.0 * (
                abs(Dumux::tripleProduct(c(3)-c(1), c(4)-c(0), c(5)-c(0)))
                + abs(Dumux::tripleProduct(c(5)-c(2), c(0)-c(2), c(1)-c(2)))
            );
        }
        else
            return std::numeric_limits<ctype>::quiet_NaN();
    }
    else
        return std::numeric_limits<ctype>::quiet_NaN();
}

/*!
 * \ingroup Geometry
 * \brief Overload for random access containers storing the corners
 */
template<int dim, std::ranges::random_access_range Corners>
auto convexPolytopeVolume(Dune::GeometryType type, const Corners& c)
{ return convexPolytopeVolume<dim>(type, [&] (unsigned int i) { return c[i]; }); }

/*!
 * \ingroup Geometry
 * \brief The volume of a given geometry
 */
template<class Geometry>
auto convexPolytopeVolume(const Geometry& geo)
{
    const auto v = convexPolytopeVolume<Geometry::mydimension>(
        geo.type(), [&](unsigned int i){ return geo.corner(i); }
    );

    // fall back to the method of the geometry if no specialized
    // volume function is implemented for the geometry type
    return std::isnan(v) ? geo.volume() : v;
}

/*!
 * \ingroup Geometry
 * \brief The volume of a given geometry
 */
template<class Geometry>
auto volume(const Geometry& geo, unsigned int integrationOrder = 4)
{
    double volume = 0.0;
    const auto rule = Dune::QuadratureRules<double, Geometry::mydimension>::rule(geo.type(), integrationOrder);
    for (const auto& qp : rule)
        volume += geo.integrationElement(qp.position())*qp.weight();
    return volume;
}

} // end namespace Dumux

#endif
