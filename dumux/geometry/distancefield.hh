// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Geometry
 * \copydoc Dumux::DistanceField
 */
#ifndef DUMUX_GEOMETRY_DISTANCE_FIELD_HH
#define DUMUX_GEOMETRY_DISTANCE_FIELD_HH

#include <memory>
#include <vector>
#include <utility>

#include <dumux/geometry/distance.hh>
#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/geometry/geometricentityset.hh>

namespace Dumux {

/*!
 * \ingroup Geometry
 * \brief Class to calculate the closest distance from a point to a given set of geometries describing the domain's boundaries.
 *        Internally uses an AABB tree representation of the geometries for logarithmic distance queries
 * \tparam Geometry The (dune) geometry type.
 */
template<class Geometry>
class AABBDistanceField
{
    using Point = typename Geometry::GlobalCoordinate;
    using Scalar = typename Geometry::ctype;
    using GeoSet = GeometriesEntitySet<Geometry>;
    using AABBTree = BoundingBoxTree<GeoSet>;

    class PointGeometry
    {
    public:
        static constexpr int mydimension = 0;
        static constexpr int coorddimension = Geometry::coorddimension;
        using ctype = typename Geometry::ctype;

        PointGeometry(Point&& point) : point_(std::move(point)) {}
        const Point& corner(std::size_t i) const { assert(i == 0); return point_; }
        std::size_t corners() const { return 1; }
    private:
        Point point_;
    };

    using PointGeoSet = GeometriesEntitySet<PointGeometry>;
    using AABBTreeMidPoints = BoundingBoxTree<PointGeoSet>;

public:
    /*!
     * \brief The constructor.
     * \param geometries A vector of geometries describing the boundaries of the spatial domain.
     */
    AABBDistanceField(const std::vector<Geometry>& geometries)
    : tree_(std::make_unique<AABBTree>(std::make_shared<GeoSet>(geometries)))
    {
        std::vector<PointGeometry> points;
        points.reserve(geometries.size());
        for (const auto& geo : geometries)
        {
            auto center = geo.center();
            points.emplace_back(std::move(center));
        }

        pointTree_ = std::make_unique<AABBTreeMidPoints>(
            std::make_shared<PointGeoSet>(std::move(points))
        );
    }

    /*!
     * \brief Returns the distance from a point to the closest geometry on the domain's boundary, as well as the index
     *        of the closest geometry.
     * \param p The location at which the closest distance is evaluated.
     */
    std::pair<Scalar, std::size_t> distanceAndIndex(const Point& p) const
    { return distanceAndIndex_(p); }

    /*!
     * \brief Returns the distance from a point to the closest geometry on the domain's boundary.
     * \param p The location at which the closest distance is evaluated.
     */
    Scalar distance(const Point& p) const
    { return distanceAndIndex_(p).first; }

private:
    std::pair<Scalar, std::size_t> distanceAndIndex_(const Point& p) const
    {
        // find a good first guess by checking the mid point tree
        const auto minSquaredDistanceEstimate = squaredDistance(p, *pointTree_);

        // find actual entity and distance to the entity's geometry
        // we choose the distance estimate a bit larger in case it actually is already the minimum distance
        const auto [squaredDistance, index] = closestEntity(p, *tree_, minSquaredDistanceEstimate*1.00001);

        using std::sqrt;
        return { sqrt(squaredDistance), index };
    }

    std::unique_ptr<AABBTree> tree_;
    std::unique_ptr<AABBTreeMidPoints> pointTree_;
};

/*!
 * \ingroup Geometry
 * \brief Class to calculate the closest distance from a point to a given set of geometries describing the domain's boundaries.
 * \tparam Geometry The (dune) geometry type.
 * \note Defaults to Dumux::AABBDistanceField
 */
template<class Geometry>
using DistanceField = AABBDistanceField<Geometry>;

} // end namespace Dumux

#endif
