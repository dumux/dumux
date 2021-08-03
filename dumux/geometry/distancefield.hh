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
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Geometry
 * \copydoc Dumux::DistanceField
 */
#ifndef DUMUX_GEOMETRY_DISTANCE_FIELD_HH
#define DUMUX_GEOMETRY_DISTANCE_FIELD_HH

#include <vector>
#include <array>
#include <dune/common/fvector.hh>
#include <dune/common/reservedvector.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dumux/common/enumerate.hh>
#include <dumux/geometry/distance.hh>
#include <dumux/geometry/circumsphere.hh>

namespace Dumux {

/*!
 * \ingroup Geometry
 * \brief Class to calculate the closest distance from a point to a given set of geometries describing the domain's boundaries.
 *        Uses a brute force search algorithm optionally enhanced by a check for bounding spheres.
 * \tparam Geometry The (dune) geometry type.
 * \tparam useBoundingSpheres Determines whether bounding spheres around the boundary geometries are constructed
 *         which may help to speed up the algorithm. If the boundary only consist of very few geometries (< 100),
 *         switching off this mechanism may actually be faster.
 */
template<class Geometry, bool useBoundingSpheres = true>
class DistanceField
{
    using Point = typename Geometry::GlobalCoordinate;
    using Scalar = typename Geometry::ctype;

    struct TriangleMLGTraits : public Dune::MultiLinearGeometryTraits<Scalar>
    {
        // we use static vectors to store the corners as we know
        // the number of corners in advance
        template<int mydim, int cdim>
        struct CornerStorage
        {
            using Type = std::array<Dune::FieldVector<Scalar, cdim>, 3>;
        };

        // we know all scvfs will have the same geometry type
        template<int mydim>
        struct hasSingleGeometryType
        {
            static const bool v = true;
            static const unsigned int topologyId = Dune::GeometryTypes::simplex(mydim).id();
        };
    };


    static constexpr auto dim = Geometry::mydimension;
    static constexpr auto dimWorld = Geometry::coorddimension;

    using Triangle = Dune::MultiLinearGeometry<Scalar, 2, dimWorld,
                                               TriangleMLGTraits>;

    using BoundingSphere = Sphere<Scalar, dimWorld>;
    using BoundingSpheresStorage = std::conditional_t<(dim == 1),
                                                      std::vector<BoundingSphere>,
                                                      std::vector<Dune::ReservedVector<BoundingSphere, 2>>>;

    static_assert(dimWorld > 1, "Distances cannot be computed for 1D domains.");

public:

    /*!
     * \brief The constructor.
     * \param geometries A vector of geometries describing the boundaries of the spatial domain.
     * \note the caller has to make sure that the object stays alive for the lifetime of the distance field
     */
    DistanceField(const std::vector<Geometry>& geometries)
    : geometries_(Dune::stackobject_to_shared_ptr(geometries))
    { initialize_(); }

    /*!
     * \brief Returns the distance from a point to the closest geometry on the domain's boundary, as well as the index
     *        of the closest geometry.
     * \param p The location at which the closest distance is evaluated.
     */
    std::pair<Scalar, std::size_t> distanceAndIndex(const Point& p) const
    {
        return getDistanceAndIndex_(p, std::integral_constant<bool, useBoundingSpheres>{});
    }

    /*!
     * \brief Returns the distance from a point to the closest geometry on the domain's boundary.
     * \param p The location at which the closest distance is evaluated.
     */
    Scalar distance(const Point& p) const
    {
        return distanceAndIndex(p).first;
    }

private:

    /*!
     * \brief Prepares the bounding spheres used for speeding up the search for the closest distance.
     */
    void initialize_()
    {
        // return early if bounding spheres are not used
        if constexpr (!useBoundingSpheres)
            return;

        boundingSpheres_.clear();
        boundingSpheres_.resize(geometries_->size());
        Point boundingSpheresCenter(0.0);
        std::size_t numSpheres = 0;

        for (auto&& [idx, geo] : enumerate(*geometries_))
        {
            if constexpr (dim == 1) // geometry is segment
            {
                const Scalar radius = (geo.corner(0) - geo.corner(1)).two_norm();
                auto center = 0.5*(geo.corner(0) + geo.corner(1));
                boundingSpheresCenter += center;
                ++numSpheres;
                boundingSpheres_[idx] = BoundingSphere{std::move(center), radius};
            }
            else // geometry is triangle or quadrilateral
            {
                static_assert(dim == 2, "Geometry must be two-dimensional");
                if (geo.corners() == 3)
                {
                    // object is already triangle
                    const auto sphere = circumSphereOfTriangle(geo);
                    boundingSpheresCenter += sphere.center();
                    ++numSpheres;
                    boundingSpheres_[idx].push_back(std::move(sphere));
                }
                else if (geo.corners() == 4)
                {
                    // object is quadrilateral, split into two triangles
                    const std::array<Triangle, 2> triangles{Triangle(Dune::GeometryTypes::simplex(2),
                                                                     std::array{geo.corner(0), geo.corner(1), geo.corner(2)}),
                                                            Triangle(Dune::GeometryTypes::simplex(2),
                                                                     std::array{geo.corner(0), geo.corner(2), geo.corner(3)})};

                    for (const auto& tria : triangles)
                    {
                        const auto sphere = circumSphereOfTriangle(tria);
                        boundingSpheresCenter += sphere.center();
                        ++numSpheres;
                        boundingSpheres_[idx].push_back(std::move(sphere));
                    }
                }
                else
                    DUNE_THROW(Dune::NotImplemented, "Object with " << geo.corners() << " not supported");
            }
        }

        boundingSpheresCenter /= numSpheres;

        Scalar squaredMinDistToCenter = std::numeric_limits<Scalar>::max();
        std::size_t closestSphereToCenterIdx = 0;

        for (auto&& [idx, geo] : enumerate(*geometries_))
        {
            if constexpr (dim == 1) // geometry is segment
            {
                if (const Scalar d = (boundingSpheres_[idx].center() - boundingSpheresCenter).two_norm2(); d < squaredMinDistToCenter)
                {
                    squaredMinDistToCenter = d;
                    closestSphereToCenterIdx = idx;
                }
            }
            else // geometry is triangle or quadrilateral
            {
                for (const auto& sphere : boundingSpheres_[idx])
                {
                    if (const Scalar d = (sphere.center() - boundingSpheresCenter).two_norm2(); d < squaredMinDistToCenter)
                    {
                        squaredMinDistToCenter = d;
                        closestSphereToCenterIdx = idx;
                    }
                }
            }
        }

        centerGeometryIndex_ = closestSphereToCenterIdx;
    }

    // overload in case bounding spheres are used
    std::pair<Scalar, std::size_t> getDistanceAndIndex_(const Point& p, std::true_type) const
    {
        using std::min;

        // Find a good initial value for the minimum distance. We chose the distance from p to the geometry
        // object closest to the mean value of the bounding sphere's centers.
        std::size_t indexOfClosestEntity = centerGeometryIndex_;
        Scalar minDistance = [&]()
        {
            const auto& centerGeometry = (*geometries_)[centerGeometryIndex_];
            if constexpr (dim == 1)
                return distancePointSegment(p, centerGeometry.corner(0), centerGeometry.corner(1));
            else
            {
                if (centerGeometry.corners() == 3)
                    return distancePointTriangle(p, centerGeometry.corner(0), centerGeometry.corner(1), centerGeometry.corner(2));
                else
                    return min(distancePointTriangle(p, centerGeometry.corner(0), centerGeometry.corner(1), centerGeometry.corner(2)),
                            distancePointTriangle(p, centerGeometry.corner(0), centerGeometry.corner(2), centerGeometry.corner(3)));
            }
        }();

        for (auto&& [idx, geo] : enumerate(*geometries_))
        {
            // We skip the check of a geometry if the closest distance from p to its circumsphere d is greater
            // than the current minimum distance:
            // 1.) d = distance(p, circumsphere) = distance(p, centerOfCircumSphere) - radiusOfCircumSphere
            // 2.) Skip if d > minDistance   --> we want avoid the rather expensive calculation of d (contains sqrt)
            // 3.) Reformulate: Skip if d^2 > minDistance^2
            // 4.) distance(p, centerOfCircumSphere)^2 - radiusOfCircumSphere^2 > minDistance^2
            // 5.) distance(p, centerOfCircumSphere)^2 > minDistance^2 + radiusOfCircumSphere^2
            const auto skipGeometry = [&](const auto& boundingSphere)
            {
                const auto squaredDistancePointSphereCenter = (p - boundingSphere.center()).two_norm2();
                const auto squaredMinDist = minDistance*minDistance;
                const auto squaredRadius = boundingSphere.radius()*boundingSphere.radius();
                return squaredDistancePointSphereCenter > squaredRadius &&
                        squaredDistancePointSphereCenter > (squaredMinDist
                                                            + minDistance*boundingSphere.radius()
                                                            + squaredRadius);
            };

            if constexpr (dim == 1)
            {
                if (!skipGeometry(boundingSpheres_[idx]))
                {
                    if (const Scalar d = distancePointSegment(p, geo.corner(0), geo.corner(1)); d < minDistance)
                    {
                        minDistance = d;
                        indexOfClosestEntity = idx;
                    }
                }
            }
            else
            {
                if (!skipGeometry(boundingSpheres_[idx][0]))
                {
                    if (const Scalar d = distancePointTriangle(p, geo.corner(0), geo.corner(1), geo.corner(2)); d < minDistance)
                    {
                        minDistance = d;
                        indexOfClosestEntity = idx;
                    }
                }

                // Check the second triangle in quadrilateral objects.
                if (geo.corners() == 4 && !skipGeometry(boundingSpheres_[idx][1]))
                {
                    if (const Scalar d = distancePointTriangle(p, geo.corner(0), geo.corner(2), geo.corner(3)); d < minDistance)
                    {
                        minDistance = d;
                        indexOfClosestEntity = idx;
                    }
                }
            }
        }

        return { minDistance, indexOfClosestEntity };
    }

    // overload in case bounding spheres are not used
    std::pair<Scalar, std::size_t> getDistanceAndIndex_(const Point& p, std::false_type) const
    {
        using std::min;
        Scalar minSquaredDistance = std::numeric_limits<Scalar>::max();
        std::size_t indexOfClosestEntity = 0;

        for (auto&& [idx, geo] : enumerate(*geometries_))
        {
            if constexpr (dim == 1)
            {
                if (const Scalar d = squaredDistancePointSegment(p, geo.corner(0), geo.corner(1)); d < minSquaredDistance)
                {
                    minSquaredDistance = d;
                    indexOfClosestEntity = idx;
                }
            }
            else
            {
                if (const Scalar d = squaredDistancePointTriangle(p, geo.corner(0), geo.corner(1), geo.corner(2)); d < minSquaredDistance)
                {
                    minSquaredDistance = d;
                    indexOfClosestEntity = idx;
                }

                if (geo.corners() == 4)
                {
                    if (const Scalar d = squaredDistancePointTriangle(p, geo.corner(0), geo.corner(2), geo.corner(3)); d < minSquaredDistance)
                    {
                        minSquaredDistance = d;
                        indexOfClosestEntity = idx;
                    }
                }
            }
        }

        using std::sqrt;
        return { sqrt(minSquaredDistance), indexOfClosestEntity };
    }

    std::shared_ptr<const std::vector<Geometry>> geometries_;
    BoundingSpheresStorage boundingSpheres_;
    std::size_t centerGeometryIndex_;
};

} // end namespace Dumux

#endif
