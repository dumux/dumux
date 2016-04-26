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
 * \brief A class for collision detection of two geometries
 *        and computation of intersection corners
 */
#ifndef DUMUX_GEOMETRY_COLLISION_HH
#define DUMUX_GEOMETRY_COLLISION_HH

#include <dune/common/exceptions.hh>
#include <dune/common/promotiontraits.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/common/math.hh>

namespace Dumux
{

//! A class for geometry collision detection and intersection calculation
template
<class Geometry1, class Geometry2,
  int dimworld = Geometry1::coorddimension,
  int dim1 = Geometry1::mydimension,
  int dim2 = Geometry2::mydimension>
class GeometryCollision
{
    using Scalar = typename Dune::PromotionTraits<typename Geometry1::ctype, typename Geometry2::ctype>::PromotedType;
    using GlobalPosition = Dune::FieldVector<Scalar, dimworld>;

public:
    //! Determine if the two geometries intersect and compute the intersection corners
    static bool collide(const Geometry1& geo1, const Geometry2& geo2, std::vector<GlobalPosition>& intersection)
    {
        static_assert(dimworld == Geometry2::coorddimension, "Can only collide geometries of same coordinate dimension");
        DUNE_THROW(Dune::NotImplemented, "Geometry collision detection with intersection computation not implemented for dimworld = "
                                          << dimworld << ", dim1 = " << dim1 << ", dim2 = " << dim2);
    }
};

//! Geometry collision detection with 3d and 1d geometry in 3d space
template <class Geometry1, class Geometry2>
class GeometryCollision<Geometry1, Geometry2, 3, 3, 1>
{
    static const int dimworld = 3;
    static const int dim1 = 3;
    static const int dim2 = 1;
    static const int dimis = dim2; // the intersection dimension
    using Scalar = typename Dune::PromotionTraits<typename Geometry1::ctype, typename Geometry2::ctype>::PromotedType;
    using ReferenceElements = typename Dune::ReferenceElements<Scalar, dim1>;
    using GlobalPosition = Dune::FieldVector<Scalar, dimworld>;

    static constexpr Scalar eps_ = 1.5e-7; // base epsilon for floating point comparisons

public:
    /*!
     *  \brief Colliding segment and convex polyhedron
     *  \note Algorithm from "Real-Time Collision Detection" by Christer Ericson
     *        Basis is the theorem that for any two non-intersecting convex polyhedrons
     *        a separating plane exists.
     *  \param intersection If the geometries collide intersection holds the corner points of
     *        the intersection object in global coordinates.
     */
    static bool collide(const Geometry1& geo1, const Geometry2& geo2, std::vector<GlobalPosition>& intersection)
    {
        static_assert(dimworld == Geometry2::coorddimension, "Can only collide geometries of same coordinate dimension");

        const auto a = geo2.corner(0);
        const auto b = geo2.corner(1);
        const auto d = b - a;

        // The initial interval is the whole segment
        // afterward we start clipping the interval
        // by the planes decribed by the facet
        Scalar tfirst = 0;
        Scalar tlast = 1;

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
                facets = {{1, 2, 0}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3}};
                break;
            default:
                DUNE_THROW(Dune::NotImplemented, "Collision of segment and geometry of type " << geo1.type() << ", "<< geo1.corners() << " corners.");
        }

        for (const auto& f : facets)
        {
            // compute normal vector by cross product
            const auto v0 = geo1.corner(f[1]) - geo1.corner(f[0]);
            const auto v1 = geo1.corner(f[2]) - geo1.corner(f[0]);
            const auto eps = eps_*v0.two_norm();

            auto n = Dumux::crossProduct(v0, v1);
            n /= n.two_norm();

            const Scalar denom = n*d;
            const Scalar dist = n*(a-geo1.corner(f[0]));

            // if denominator is zero the segment in parallel to
            // the plane. If the distance is positive there is no intersection
            if (std::abs(denom) < eps)
            {
                if (dist > eps)
                    return false;
            }
            else // not parallel: compute line-plane intersection
            {
                const Scalar t = -dist / denom;
                // if entering half space cut tfirst if t is larger
                if (std::signbit(denom))
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
        intersection = {geo2.global(tfirst), geo2.global(tlast)};
        return true;
    }
};

} // end namespace Dumux

# endif
