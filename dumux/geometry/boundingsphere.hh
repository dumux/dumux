// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Geometry
 * \brief A function to compute bounding spheres of points clouds or convex polytopes
 */
#ifndef DUMUX_GEOMETRY_BOUNDINGSPHERE_HH
#define DUMUX_GEOMETRY_BOUNDINGSPHERE_HH

#include <algorithm>
#include <vector>

#include <dumux/geometry/sphere.hh>

namespace Dumux {

/*!
 * \ingroup Geometry
 * \brief Computes a bounding sphere of a convex polytope geometry (Dune::Geometry interface)
 * \note The bounding sphere is not necessarily the minimum enclosing sphere
 *       but computation is fast (AABB method)
 */
template<class ConvexGeometry>
static inline Sphere<typename ConvexGeometry::ctype, ConvexGeometry::coorddimension>
boundingSphere(const ConvexGeometry& geometry)
{
    constexpr int dimWorld = ConvexGeometry::coorddimension;
    assert(geometry.corners() >= 1);

    auto corner = geometry.corner(0);
    auto xMin = corner, xMax = corner;
    using std::max; using std::min;

    // Compute the min and max over the remaining vertices
    for (std::size_t i = 1; i < geometry.corners(); ++i)
    {
        corner = geometry.corner(i);
        for (std::size_t dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
        {
            xMin[dimIdx] = min(xMin[dimIdx], corner[dimIdx]);
            xMax[dimIdx] = max(xMax[dimIdx], corner[dimIdx]);
        }
    }

    auto center = 0.5*xMax + xMin;
    auto radius = (center-xMax).two_norm();
    return { std::move(center), std::move(radius) };
}

namespace Detail {
template<class Scalar, int dim>
struct PointsToGeometryWrapper
{
    const std::vector<Dune::FieldVector<Scalar, dim>>& points_;
    PointsToGeometryWrapper(const std::vector<Dune::FieldVector<Scalar, dim>>& points)
    : points_(points) {}

    static constexpr int coorddimension = dim;
    using ctype = Scalar;

    auto corners() const { return points_.size(); }
    const auto& corner(std::size_t i) const { return points_[i]; }
};
} // end namespace Detail

/*!
 * \ingroup Geometry
 * \brief Computes a bounding sphere of points
 * \note The bounding sphere is not necessarily the minimum enclosing sphere
 *       but computation is fast (AABB method)
 */
template<class Scalar, int dim>
static inline Sphere<Scalar, dim>
boundingSphere(const std::vector<Dune::FieldVector<Scalar, dim>>& points)
{
    Detail::PointsToGeometryWrapper<Scalar, dim> geometry(points);
    return boundingSphere(geometry);
}

} // end namespace Dumux

#endif
