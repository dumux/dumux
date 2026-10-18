// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Contact
 * \brief Closest point on a segment or triangle with its barycentric weights
 */
#ifndef DUMUX_SOLIDMECHANICS_CONTACT_CLOSEST_POINT_HH
#define DUMUX_SOLIDMECHANICS_CONTACT_CLOSEST_POINT_HH

#include <algorithm>
#include <array>

#include <dune/common/fvector.hh>

namespace Dumux::Contact {

/*!
 * \ingroup Contact
 * \brief The closest point \f$ \mathbf{y} = \sum_j w_j \mathbf{c}_j \f$ on a simplex with corners
 *        \f$ \mathbf{c}_j \f$
 *
 * The weights of the corners not spanning the closest feature (vertex, edge or the simplex
 * itself) are exactly zero.
 */
template<class Scalar, int dimWorld, std::size_t numCorners>
struct ClosestPoint
{
    Dune::FieldVector<Scalar, dimWorld> point;
    std::array<Scalar, numCorners> weights;
    Scalar squaredDistance;
};

/*!
 * \ingroup Contact
 * \brief The closest point on the segment with the corners a and b
 */
template<class Scalar, int dimWorld>
ClosestPoint<Scalar, dimWorld, 2>
closestPointOnSegment(const Dune::FieldVector<Scalar, dimWorld>& p,
                      const Dune::FieldVector<Scalar, dimWorld>& a,
                      const Dune::FieldVector<Scalar, dimWorld>& b)
{
    const auto ab = b - a;
    const Scalar t = std::clamp((p - a)*ab/ab.two_norm2(), Scalar(0.0), Scalar(1.0));
    ClosestPoint<Scalar, dimWorld, 2> result{a, {1.0 - t, t}, 0.0};
    result.point.axpy(t, ab);
    result.squaredDistance = (p - result.point).two_norm2();
    return result;
}

/*!
 * \ingroup Contact
 * \brief The closest point on the triangle with the corners a, b and c
 * \note Region classification by Voronoi regions of the features, see Ericson,
 *       Real-Time Collision Detection (2005), section 5.1.5
 */
template<class Scalar, int dimWorld>
ClosestPoint<Scalar, dimWorld, 3>
closestPointOnTriangle(const Dune::FieldVector<Scalar, dimWorld>& p,
                       const Dune::FieldVector<Scalar, dimWorld>& a,
                       const Dune::FieldVector<Scalar, dimWorld>& b,
                       const Dune::FieldVector<Scalar, dimWorld>& c)
{
    const auto weights = [&]() -> std::array<Scalar, 3> {
        const auto ab = b - a;
        const auto ac = c - a;
        const auto ap = p - a;
        const Scalar d1 = ab*ap, d2 = ac*ap;
        if (d1 <= 0.0 && d2 <= 0.0)
            return {1.0, 0.0, 0.0};

        const auto bp = p - b;
        const Scalar d3 = ab*bp, d4 = ac*bp;
        if (d3 >= 0.0 && d4 <= d3)
            return {0.0, 1.0, 0.0};

        const Scalar vc = d1*d4 - d3*d2;
        if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
        {
            const Scalar v = d1/(d1 - d3);
            return {1.0 - v, v, 0.0};
        }

        const auto cp = p - c;
        const Scalar d5 = ab*cp, d6 = ac*cp;
        if (d6 >= 0.0 && d5 <= d6)
            return {0.0, 0.0, 1.0};

        const Scalar vb = d5*d2 - d1*d6;
        if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
        {
            const Scalar w = d2/(d2 - d6);
            return {1.0 - w, 0.0, w};
        }

        const Scalar va = d3*d6 - d5*d4;
        if (va <= 0.0 && d4 - d3 >= 0.0 && d5 - d6 >= 0.0)
        {
            const Scalar w = (d4 - d3)/((d4 - d3) + (d5 - d6));
            return {0.0, 1.0 - w, w};
        }

        const Scalar denominator = 1.0/(va + vb + vc);
        const Scalar v = vb*denominator, w = vc*denominator;
        return {1.0 - v - w, v, w};
    }();

    ClosestPoint<Scalar, dimWorld, 3> result{a, weights, 0.0};
    result.point *= weights[0];
    result.point.axpy(weights[1], b);
    result.point.axpy(weights[2], c);
    result.squaredDistance = (p - result.point).two_norm2();
    return result;
}

/*!
 * \ingroup Contact
 * \brief The closest point on a segment (two corners) or a triangle (three corners)
 */
template<class Scalar, int dimWorld, std::size_t numCorners>
ClosestPoint<Scalar, dimWorld, numCorners>
closestPointOnSimplex(const Dune::FieldVector<Scalar, dimWorld>& p,
                      const std::array<Dune::FieldVector<Scalar, dimWorld>, numCorners>& corners)
{
    static_assert(numCorners == 2 || numCorners == 3, "Only segments and triangles are implemented");
    if constexpr (numCorners == 2)
        return closestPointOnSegment(p, corners[0], corners[1]);
    else
        return closestPointOnTriangle(p, corners[0], corners[1], corners[2]);
}

} // end namespace Dumux::Contact

#endif
