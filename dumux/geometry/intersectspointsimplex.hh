// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Geometry
 * \brief Detect if a point intersects a simplex (including boundary)
 */
#ifndef DUMUX_GEOMETRY_INTERSECTS_POINT_SIMPLEX_HH
#define DUMUX_GEOMETRY_INTERSECTS_POINT_SIMPLEX_HH

#include <cmath>
#include <dune/common/fvector.hh>
#include <dumux/common/math.hh>

namespace Dumux {

/*!
 * \ingroup Geometry
 * \brief Find out whether a point is inside the tetrahedron (p0, p1, p2, p3) (dimworld is 3)
 */
template<class ctype, int dimworld, typename std::enable_if_t<(dimworld == 3), int> = 0>
bool intersectsPointSimplex(const Dune::FieldVector<ctype, dimworld>& point,
                            const Dune::FieldVector<ctype, dimworld>& p0,
                            const Dune::FieldVector<ctype, dimworld>& p1,
                            const Dune::FieldVector<ctype, dimworld>& p2,
                            const Dune::FieldVector<ctype, dimworld>& p3)
{
    // Algorithm from http://www.blackpawn.com/texts/pointinpoly/
    // See also "Real-Time Collision Detection" by Christer Ericson.
    using GlobalPosition = Dune::FieldVector<ctype, dimworld>;
    static constexpr ctype eps_ = 1.0e-7;

    // put the tetrahedron points in an array
    const GlobalPosition *p[4] = {&p0, &p1, &p2, &p3};

    // iterate over all faces
    for (int i = 0; i < 4; ++i)
    {
        // compute all the vectors from vertex (local index 0) to the other points
        const GlobalPosition v1 = *p[(i + 1)%4] - *p[i];
        const GlobalPosition v2 = *p[(i + 2)%4] - *p[i];
        const GlobalPosition v3 = *p[(i + 3)%4] - *p[i];
        const GlobalPosition v = point - *p[i];
        // compute the normal to the facet (cross product)
        GlobalPosition n1 = crossProduct(v1, v2);
        n1 /= n1.two_norm();
        // find out on which side of the plane v and v3 are
        const auto t1 = n1.dot(v);
        const auto t2 = n1.dot(v3);
        // If the point is not exactly on the plane the
        // points have to be on the same side
        const auto eps = eps_ * v1.two_norm();
        if ((t1 > eps || t1 < -eps) && std::signbit(t1) != std::signbit(t2))
            return false;
    }
    return true;
}

/*!
 * \ingroup Geometry
 * \brief Find out whether a point is inside the triangle (p0, p1, p2) (dimworld is 3)
 */
template<class ctype, int dimworld, typename std::enable_if_t<(dimworld == 3), int> = 0>
bool intersectsPointSimplex(const Dune::FieldVector<ctype, dimworld>& point,
                            const Dune::FieldVector<ctype, dimworld>& p0,
                            const Dune::FieldVector<ctype, dimworld>& p1,
                            const Dune::FieldVector<ctype, dimworld>& p2)
{
    // adapted from the algorithm from from "Real-Time Collision Detection" by Christer Ericson,
    // published by Morgan Kaufmann Publishers, (c) 2005 Elsevier Inc. (Chapter 5.4.2)
    constexpr ctype eps_ = 1.0e-7;

    // compute the normal of the triangle
    const auto v1 = p0 - p2;
    auto n = crossProduct(v1, p1 - p0);
    const ctype nnorm = n.two_norm();
    const ctype eps4 = eps_*nnorm*nnorm; // compute an epsilon for later
    n /= nnorm; // normalize

    // first check if we are in the plane of the triangle
    // if not we can return early
    using std::abs;
    auto x = p0 - point;
    x /= x.two_norm(); // normalize

    if (abs(x*n) > eps_)
        return false;

    // translate the triangle so that 'point' is the origin
    const auto a = p0 - point;
    const auto b = p1 - point;
    const auto c = p2 - point;

    // compute the normal vectors for triangles P->A->B and P->B->C
    const auto u = crossProduct(b, c);
    const auto v = crossProduct(c, a);

    // they have to point in the same direction or be orthogonal
    if (u*v < 0.0 - eps4)
        return false;

    // compute the normal vector for triangle P->C->A
    const auto w = crossProduct(a, b);

    // it also has to point in the same direction or be orthogonal
    if (u*w < 0.0 - eps4)
        return false;

    // check if point is on the line of one of the edges
    if (u.two_norm2() < eps4)
        return b*c < 0.0 + eps_*nnorm;
    if (v.two_norm2() < eps4)
        return a*c < 0.0 + eps_*nnorm;

    // now the point must be in the triangle (or on the faces)
    return true;
}

/*!
 * \ingroup Geometry
 * \brief Find out whether a point is inside the triangle (p0, p1, p2) (dimworld is 2)
 */
template<class ctype, int dimworld, typename std::enable_if_t<(dimworld == 2), int> = 0>
bool intersectsPointSimplex(const Dune::FieldVector<ctype, dimworld>& point,
                            const Dune::FieldVector<ctype, dimworld>& p0,
                            const Dune::FieldVector<ctype, dimworld>& p1,
                            const Dune::FieldVector<ctype, dimworld>& p2)
{
    static constexpr ctype eps_ = 1.0e-7;

    // Use barycentric coordinates
    const ctype A = 0.5*(-p1[1]*p2[0] + p0[1]*(p2[0] - p1[0])
                         +p1[0]*p2[1] + p0[0]*(p1[1] - p2[1]));
    const ctype sign = std::copysign(1.0, A);
    const ctype s = sign*(p0[1]*p2[0] + point[0]*(p2[1]-p0[1])
                         -p0[0]*p2[1] + point[1]*(p0[0]-p2[0]));
    const ctype t = sign*(p0[0]*p1[1] + point[0]*(p0[1]-p1[1])
                         -p0[1]*p1[0] + point[1]*(p1[0]-p0[0]));
    const ctype eps = sign*A*eps_;

    return (s > -eps
            && t > -eps
            && (s + t) < 2*A*sign + eps);
}

/*!
 * \ingroup Geometry
 * \brief Find out whether a point is inside the interval (p0, p1) (dimworld is 2 or 3)
 * \note We assume the given interval has non-zero length and use it to scale the epsilon
 */
template<class ctype, int dimworld, typename std::enable_if_t<(dimworld == 3 || dimworld == 2), int> = 0>
bool intersectsPointSimplex(const Dune::FieldVector<ctype, dimworld>& point,
                            const Dune::FieldVector<ctype, dimworld>& p0,
                            const Dune::FieldVector<ctype, dimworld>& p1)
{
    using GlobalPosition = Dune::FieldVector<ctype, dimworld>;
    static constexpr ctype eps_ = 1.0e-7;

    // compute the vectors between p0 and the other points
    const GlobalPosition v1 = p1 - p0;
    const GlobalPosition v2 = point - p0;

    const ctype v1norm = v1.two_norm();
    const ctype v2norm = v2.two_norm();

    // early exit if point and p0 are the same
    if (v2norm < v1norm*eps_)
        return true;

    // early exit if the point is outside the segment (no epsilon in the
    // first statement because we already did the above equality check)
    if (v1.dot(v2) < 0.0 || v2norm > v1norm*(1.0 + eps_))
        return false;

    // If the area spanned by the 2 vectors is zero, the points are colinear.
    // If that is the case, the given point is on the segment.
    const auto n = crossProduct(v1, v2);
    const auto eps2 = v1norm*v1norm*eps_;
    if constexpr (dimworld == 3)
        return n.two_norm2() < eps2*eps2;
    else
    {
        using std::abs;
        return abs(n) < eps2;
    }
}

/*!
 * \ingroup Geometry
 * \brief Find out whether a point is inside the interval (p0, p1) (dimworld is 1)
 */
template<class ctype, int dimworld, typename std::enable_if_t<(dimworld == 1), int> = 0>
bool intersectsPointSimplex(const Dune::FieldVector<ctype, dimworld>& point,
                            const Dune::FieldVector<ctype, dimworld>& p0,
                            const Dune::FieldVector<ctype, dimworld>& p1)
{
    static constexpr ctype eps_ = 1.0e-7;

    // sort the interval so interval[1] is the end and interval[0] the start
    const ctype *interval[2] = {&p0[0], &p1[0]};
    if (*interval[0] > *interval[1])
        std::swap(interval[0], interval[1]);

    const ctype v1 = point[0] - *interval[0];
    const ctype v2 = *interval[1] - *interval[0]; // always positive

    // the coordinates are the same
    using std::abs;
    if (abs(v1) < v2*eps_)
        return true;

    // the point doesn't coincide with p0
    // so if p0 and p1 are equal it's not inside
    if (v2 < 1.0e-30)
        return false;

    // the point is inside if the length is
    // smaller than the interval length and the
    // sign of v1 & v2 are the same
    using std::signbit;
    return (!signbit(v1) && abs(v1) < v2*(1.0 + eps_));
}

} // end namespace Dumux

#endif
