// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkModels
 * \brief Geometry of the pores and throats of a regular triangulation of a sphere packing
 *
 * A pore is a tetrahedron of the regular (weighted Delaunay) triangulation of the sphere centres
 * with weights r^2, a throat is a triangular facet shared by two tetrahedra. The definitions follow
 * the pore-scale finite volume method \cite Chareyre2012. The solid part of a region is approximated
 * by the sphere sectors inside it, i.e. spheres reaching beyond the opposite face of a tetrahedron are
 * not cut.
 */
#ifndef DUMUX_PNM_EXTRACTION_SPHERE_PACKING_GEOMETRY_HH
#define DUMUX_PNM_EXTRACTION_SPHERE_PACKING_GEOMETRY_HH

#include <algorithm>
#include <array>
#include <cmath>
#include <utility>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Dumux::PoreNetwork::SpherePacking {

template<class Scalar>
using Point = Dune::FieldVector<Scalar, 3>;

namespace Detail {

template<class Scalar>
Point<Scalar> cross(const Point<Scalar>& a, const Point<Scalar>& b)
{ return {a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]}; }

template<class Scalar>
Scalar tripleProduct(const Point<Scalar>& a, const Point<Scalar>& b, const Point<Scalar>& c)
{ return a*cross(b, c); }

//! Smallest positive root of A x^2 + 2 B x + C = 0, else the largest real root, else zero
template<class Scalar>
Scalar selectRoot(Scalar A, Scalar B, Scalar C)
{
    using std::abs; using std::sqrt; using std::max; using std::min;
    if (abs(A) < 1e-12)
        return B != 0.0 ? -C/(2.0*B) : 0.0;

    const Scalar discriminant = B*B - A*C;
    if (discriminant < 0.0)
        return 0.0;

    // stable form avoiding cancellation between -B and the square root
    const Scalar q = -(B + std::copysign(sqrt(discriminant), B));
    const Scalar x1 = q/A;
    const Scalar x2 = q != 0.0 ? C/q : x1;
    const Scalar lo = min(x1, x2), hi = max(x1, x2);
    return lo > 0.0 ? lo : hi;
}

/*!
 * \brief Centre (relative to x0) and radius of the circle or sphere touching the given spheres from outside
 *
 * With the centre c = x0 + a + b rho, the tangency conditions |c - x_i| = r_i + rho reduce to
 * (c - x0).d_i = e_i + f_i rho with d_i = x_i - x0, e_i = (|d_i|^2 - r_i^2 + r0^2)/2 and f_i = r0 - r_i,
 * where a and b are the solutions for the right-hand sides e and f within the span of the d_i,
 * and the remaining condition |a + b rho| = r0 + rho is a quadratic equation for rho.
 * The vector a is the power centre relative to x0.
 */
template<class Scalar, int n>
std::pair<Point<Scalar>, Scalar> outerTangentBall(const std::array<Point<Scalar>, n>& d,
                                                  const Dune::FieldVector<Scalar, n>& e,
                                                  const Dune::FieldVector<Scalar, n>& f,
                                                  Scalar r0)
{
    Dune::FieldMatrix<Scalar, n, n> gram;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            gram[i][j] = d[i]*d[j];

    Dune::FieldVector<Scalar, n> alpha, beta;
    gram.solve(alpha, e);
    gram.solve(beta, f);

    Point<Scalar> a(0.0), b(0.0);
    for (int i = 0; i < n; ++i)
    {
        a.axpy(alpha[i], d[i]);
        b.axpy(beta[i], d[i]);
    }
    const Scalar rho = selectRoot(b.two_norm2() - 1.0, a*b - r0, a.two_norm2() - r0*r0);
    return {a.axpy(rho, b), rho};
}

} // end namespace Detail

/*!
 * \brief Solid angle of the triangle (p1, p2, p3) seen from the apex
 * \note Van Oosterom and Strackee (1983), IEEE Trans. Biomed. Eng. 30, 125-126,
 *       with atan2 so that solid angles larger than pi are correct
 */
template<class Scalar>
Scalar solidAngle(const Point<Scalar>& apex, const Point<Scalar>& p1, const Point<Scalar>& p2, const Point<Scalar>& p3)
{
    using std::abs; using std::atan2;
    const auto a = p1 - apex, b = p2 - apex, c = p3 - apex;
    const Scalar la = a.two_norm(), lb = b.two_norm(), lc = c.two_norm();
    const Scalar numerator = abs(Detail::tripleProduct(a, b, c));
    const Scalar denominator = la*lb*lc + (a*b)*lc + (a*c)*lb + (b*c)*la;
    return 2.0*atan2(numerator, denominator);
}

//! Interior angle of the triangle (apex, p1, p2) at the apex
template<class Scalar>
Scalar planeAngle(const Point<Scalar>& apex, const Point<Scalar>& p1, const Point<Scalar>& p2)
{
    using std::atan2;
    const auto a = p1 - apex, b = p2 - apex;
    return atan2(Detail::cross(a, b).two_norm(), a*b);
}

template<class Scalar>
Scalar tetrahedronVolume(const std::array<Point<Scalar>, 4>& x)
{
    using std::abs;
    return abs(Detail::tripleProduct(x[1] - x[0], x[2] - x[0], x[3] - x[0]))/6.0;
}

template<class Scalar>
Scalar triangleArea(const Point<Scalar>& x0, const Point<Scalar>& x1, const Point<Scalar>& x2)
{ return 0.5*Detail::cross(x1 - x0, x2 - x0).two_norm(); }

/*!
 * \brief Power centre of four spheres: the point c with equal power |c - x_i|^2 - r_i^2
 *
 * It is the weighted circumcentre of the tetrahedron and the vertex of the power diagram dual to it.
 */
template<class Scalar>
Point<Scalar> powerCenter(const std::array<Point<Scalar>, 4>& x, const std::array<Scalar, 4>& r)
{
    Dune::FieldMatrix<Scalar, 3, 3> m;
    Dune::FieldVector<Scalar, 3> rhs;
    for (int i = 0; i < 3; ++i)
    {
        const auto d = x[i+1] - x[0];
        m[i] = d;
        rhs[i] = 0.5*(d.two_norm2() - r[i+1]*r[i+1] + r[0]*r[0]);
    }
    Point<Scalar> c;
    m.solve(c, rhs);
    return c += x[0];
}

//! A sphere or, within a plane, a circle
template<class Scalar>
struct Ball
{
    Point<Scalar> center;
    Scalar radius;
};

/*!
 * \brief Largest sphere in the void between four spheres (3D Apollonius problem)
 * \return the sphere touching all four from outside; its radius is zero or negative if the
 *         spheres close the void
 */
template<class Scalar>
Ball<Scalar> inscribedSphere(const std::array<Point<Scalar>, 4>& x, const std::array<Scalar, 4>& r)
{
    std::array<Point<Scalar>, 3> d;
    Dune::FieldVector<Scalar, 3> e, f;
    for (int i = 0; i < 3; ++i)
    {
        d[i] = x[i+1] - x[0];
        e[i] = 0.5*(d[i].two_norm2() - r[i+1]*r[i+1] + r[0]*r[0]);
        f[i] = r[0] - r[i+1];
    }
    auto [center, radius] = Detail::outerTangentBall<Scalar, 3>(d, e, f, r[0]);
    return {center += x[0], radius};
}

/*!
 * \brief Radius of the largest sphere around the given centre that does not intersect the spheres
 *
 * Used as pore radius where no sphere touches all spheres of a flat tetrahedron from outside.
 */
template<class Scalar, class Centers, class Radii>
Scalar clearance(const Point<Scalar>& center, const Centers& x, const Radii& r)
{
    using std::min;
    Scalar radius = (center - x[0]).two_norm() - r[0];
    for (std::size_t i = 1; i < x.size(); ++i)
        radius = min(radius, (center - x[i]).two_norm() - r[i]);
    return radius;
}

/*!
 * \brief Largest circle in the plane of three sphere centres between the cross-sections of the
 *        spheres (2D Apollonius problem)
 */
template<class Scalar>
Ball<Scalar> inscribedCircle(const std::array<Point<Scalar>, 3>& x, const std::array<Scalar, 3>& r)
{
    std::array<Point<Scalar>, 2> d;
    Dune::FieldVector<Scalar, 2> e, f;
    for (int i = 0; i < 2; ++i)
    {
        d[i] = x[i+1] - x[0];
        e[i] = 0.5*(d[i].two_norm2() - r[i+1]*r[i+1] + r[0]*r[0]);
        f[i] = r[0] - r[i+1];
    }
    auto [center, radius] = Detail::outerTangentBall<Scalar, 2>(d, e, f, r[0]);
    return {center += x[0], radius};
}

//! Volume of the sphere sectors inside the tetrahedron, one per vertex
template<class Scalar>
Scalar tetrahedronSolidVolume(const std::array<Point<Scalar>, 4>& x, const std::array<Scalar, 4>& r)
{
    Scalar solid = 0.0;
    for (int i = 0; i < 4; ++i)
        solid += solidAngle(x[i], x[(i+1)%4], x[(i+2)%4], x[(i+3)%4])*r[i]*r[i]*r[i]/3.0;
    return solid;
}

//! Void volume of a tetrahedron: its volume minus the sphere sectors at its vertices
template<class Scalar>
Scalar tetrahedronVoidVolume(const std::array<Point<Scalar>, 4>& x, const std::array<Scalar, 4>& r)
{ return tetrahedronVolume(x) - tetrahedronSolidVolume(x, r); }

/*!
 * \brief Area of the circular segment of the circle (x0, r0) beyond the line through x1 and x2,
 *        zero if the foot of x0 on the line is outside the edge or the circle does not cross the line
 */
template<class Scalar>
Scalar circularSegmentBeyondEdge(const Point<Scalar>& x0, Scalar r0, const Point<Scalar>& x1, const Point<Scalar>& x2)
{
    using std::sqrt; using std::acos;
    const auto edge = x2 - x1;
    const Scalar edgeLength2 = edge.two_norm2();
    const Scalar projection = (x0 - x1)*edge;
    if (projection < 0.0 || projection > edgeLength2)
        return 0.0;

    const Scalar distance2 = Detail::cross(x0 - x1, edge).two_norm2()/edgeLength2;
    if (distance2 >= r0*r0)
        return 0.0;

    const Scalar chord = 2.0*sqrt(r0*r0 - distance2);
    const Scalar angle = 2.0*acos(sqrt(distance2)/r0);
    return 0.5*(angle*r0*r0 - chord*sqrt(distance2));
}

/*!
 * \brief Fluid area of a facet: the triangle area minus the circle sectors at its vertices
 *
 * A circle that crosses the opposite edge would be counted beyond the triangle; the circular
 * segment beyond the edge is added back for the first such circle.
 */
template<class Scalar>
Scalar facetFluidArea(const std::array<Point<Scalar>, 3>& x, const std::array<Scalar, 3>& r)
{
    using std::abs;
    Scalar area = triangleArea(x[0], x[1], x[2]);
    for (int i = 0; i < 3; ++i)
        area -= 0.5*r[i]*r[i]*planeAngle(x[i], x[(i+1)%3], x[(i+2)%3]);

    for (int i = 0; i < 3; ++i)
        if (const Scalar segment = circularSegmentBeyondEdge(x[i], r[i], x[(i+1)%3], x[(i+2)%3]); segment > 0.0)
            return abs(area + segment);

    return abs(area);
}

//! Volume, void volume and wetted solid surface of the region of a throat
template<class Scalar>
struct ThroatRegion
{
    Scalar volume = 0.0;
    Scalar voidVolume = 0.0;
    Scalar solidSurface = 0.0;
    std::array<Scalar, 3> sphereSolidSurface = {}; //!< solid surface of each of the three spheres

    //! hydraulic radius as void volume over wetted solid surface
    Scalar hydraulicRadius() const
    { return solidSurface > 0.0 ? voidVolume/solidSurface : 0.0; }
};

/*!
 * \brief The region between a facet and the power centres p1, p2 of the two adjacent tetrahedra
 *
 * The region is the union of the pyramids with the facet as base and p1, p2 as apexes (Sect. 2.2.2 in
 * \cite Chareyre2012). The line through p1 and p2 is normal to the facet. A power centre may lie
 * outside its tetrahedron, beyond the facet; its pyramid then counts negatively, so that the regions
 * of all facets of a triangulation partition its void volume and wetted solid surface.
 */
template<class Scalar>
ThroatRegion<Scalar> throatRegion(const std::array<Point<Scalar>, 3>& x, const std::array<Scalar, 3>& r,
                                  const Point<Scalar>& p1, const Point<Scalar>& p2)
{
    auto normal = Detail::cross(x[1] - x[0], x[2] - x[0]);
    if (normal*(p1 - p2) < 0.0)
        normal *= -1.0;

    ThroatRegion<Scalar> region;
    Scalar solidVolume = 0.0;
    const auto addPyramid = [&](const Point<Scalar>& apex, Scalar sign) {
        region.volume += sign*tetrahedronVolume(std::array<Point<Scalar>, 4>{x[0], x[1], x[2], apex});
        for (int i = 0; i < 3; ++i)
        {
            const Scalar omega = sign*solidAngle(x[i], x[(i+1)%3], x[(i+2)%3], apex);
            solidVolume += omega*r[i]*r[i]*r[i]/3.0;
            region.sphereSolidSurface[i] += omega*r[i]*r[i];
            region.solidSurface += omega*r[i]*r[i];
        }
    };
    addPyramid(p1, normal*(p1 - x[0]) >= 0.0 ? 1.0 : -1.0);
    addPyramid(p2, normal*(p2 - x[0]) <= 0.0 ? 1.0 : -1.0);
    region.voidVolume = region.volume - solidVolume;
    return region;
}

} // end namespace Dumux::PoreNetwork::SpherePacking

#endif
