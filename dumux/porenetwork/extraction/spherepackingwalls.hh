// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkModels
 * \brief Geometry of the pores and throats of a sphere packing next to the walls of a box
 *
 * Each of the six walls is a vertex of the regular triangulation: a sphere of very large radius whose
 * surface is the wall plane (Chareyre et al. 2012; YADE FlowEngine). Tetrahedra and facets with such
 * vertices are evaluated on the wall plane, with the projections of the sphere centres onto it.
 * The walls are slip walls: they do not count as wetted surface, and the hydraulic radius of a facet
 * touching one or two walls is multiplied by 2^(-1/4) or 2^(-1/2) as in YADE.
 */
#ifndef DUMUX_PNM_EXTRACTION_SPHERE_PACKING_WALLS_HH
#define DUMUX_PNM_EXTRACTION_SPHERE_PACKING_WALLS_HH

#include <array>
#include <cmath>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dumux/porenetwork/extraction/spherepackinggeometry.hh>

namespace Dumux::PoreNetwork::SpherePacking {

/*!
 * \brief The six walls of an axis-aligned box as vertices numSpheres + k of the triangulation
 *
 * Wall k = 2*axis + side lies at lower[axis] (side 0) or upper[axis] (side 1). Its vertex is a sphere
 * of radius farFactor times the box extent in y, centred on the outer side of the wall at the centre
 * of the face, so that its surface is the wall plane.
 */
template<class Scalar>
struct Walls
{
    Point<Scalar> lower;
    Point<Scalar> upper;
    int numSpheres = 0;
    Scalar farFactor = 5e4;

    bool isWall(int vertex) const { return vertex >= numSpheres; }
    int wall(int vertex) const { return vertex - numSpheres; }
    static int axis(int wall) { return wall/2; }
    Scalar position(int wall) const { return wall % 2 == 0 ? lower[wall/2] : upper[wall/2]; }

    Point<Scalar> inwardNormal(int wall) const
    {
        Point<Scalar> n(0.0);
        n[wall/2] = wall % 2 == 0 ? 1.0 : -1.0;
        return n;
    }

    Scalar radius() const { return farFactor*(upper[1] - lower[1]); }

    Point<Scalar> center(int wall) const
    {
        Point<Scalar> c = 0.5*(lower + upper);
        c[wall/2] = position(wall);
        return c.axpy(-radius(), inwardNormal(wall));
    }

    Point<Scalar> project(Point<Scalar> x, int wall) const
    {
        x[wall/2] = position(wall);
        return x;
    }
};

//! Centre and radius of a sphere (or of a wall vertex) of the triangulation
template<class Scalar>
struct Vertex
{
    Point<Scalar> x;
    Scalar r = 0.0;
    int wall = -1; //!< index of the wall, or -1 for a sphere
};

/*!
 * \brief Largest ball touching the given spheres from outside and lying on the inner side of the given
 *        walls, with its centre in the affine space origin + span(basis)
 *
 * The tangency conditions to the spheres, minus the one to the first sphere, and the distances to the
 * walls are linear in the centre and the radius; the remaining condition is quadratic in the radius.
 * At least one sphere must be given.
 */
template<class Scalar, int dim>
Ball<Scalar> tangentBall(const std::vector<Vertex<Scalar>>& vertices, const Walls<Scalar>& walls,
                         const Point<Scalar>& origin, const std::array<Point<Scalar>, dim>& basis)
{
    const Vertex<Scalar>* first = nullptr;
    for (const auto& v : vertices)
        if (v.wall < 0) { first = &v; break; }

    // rows: coefficients of y (dim) and rho, right-hand side
    Dune::FieldMatrix<Scalar, dim, dim> matrix;
    Dune::FieldVector<Scalar, dim> constant, slope;
    int row = 0;
    for (const auto& v : vertices)
    {
        if (&v == first)
            continue;
        if (v.wall >= 0)
        {
            const auto n = walls.inwardNormal(v.wall);
            for (int k = 0; k < dim; ++k)
                matrix[row][k] = n*basis[k];
            constant[row] = walls.position(v.wall)*n[Walls<Scalar>::axis(v.wall)] - n*origin;
            slope[row] = 1.0;
        }
        else
        {
            const auto d = v.x - first->x;
            for (int k = 0; k < dim; ++k)
                matrix[row][k] = -2.0*(d*basis[k]);
            constant[row] = v.r*v.r - first->r*first->r - d*(v.x + first->x) + 2.0*(d*origin);
            slope[row] = 2.0*(v.r - first->r);
        }
        ++row;
    }

    Dune::FieldVector<Scalar, dim> a, b;
    matrix.solve(a, constant);
    matrix.solve(b, slope);
    Point<Scalar> ca = origin, cb(0.0);
    for (int k = 0; k < dim; ++k)
    {
        ca.axpy(a[k], basis[k]);
        cb.axpy(b[k], basis[k]);
    }
    const auto e = ca - first->x;
    const Scalar rho = Detail::selectRoot(cb.two_norm2() - 1.0, e*cb - first->r, e.two_norm2() - first->r*first->r);
    return {ca.axpy(rho, cb), rho};
}

/*!
 * \brief Volume of a tetrahedron with wall vertices: the prism between the face of its spheres and the
 *        wall for one wall, two pyramids for two walls, the box corner for three walls (YADE FlowEngine)
 */
template<class Scalar>
Scalar cellBulkVolume(const std::array<Vertex<Scalar>, 4>& v, const Walls<Scalar>& walls)
{
    using std::abs;
    std::vector<Point<Scalar>> spheres;
    std::vector<int> wallIds;
    for (const auto& vertex : v)
        if (vertex.wall < 0)
            spheres.push_back(vertex.x);
        else
            wallIds.push_back(vertex.wall);

    if (wallIds.empty())
        return tetrahedronVolume(std::array<Point<Scalar>, 4>{v[0].x, v[1].x, v[2].x, v[3].x});

    if (wallIds.size() == 1)
    {
        const int c = Walls<Scalar>::axis(wallIds[0]);
        const auto normal = Detail::cross(spheres[0] - spheres[1], spheres[0] - spheres[2]);
        const Scalar mean = (spheres[0][c] + spheres[1][c] + spheres[2][c])/3.0;
        return abs(0.5*normal[c]*(mean - walls.position(wallIds[0])));
    }

    if (wallIds.size() == 2)
    {
        const int c0 = Walls<Scalar>::axis(wallIds[0]), c1 = Walls<Scalar>::axis(wallIds[1]);
        const Scalar w0 = walls.position(wallIds[0]), w1 = walls.position(wallIds[1]);
        const auto& a = spheres[0];
        const auto& b = spheres[1];
        auto as = a, bs = b;
        as[c0] = bs[c0] = w0;
        const Scalar vol1 = 0.5*Detail::cross(a - bs, b - bs)[c1]*((2.0*b[c1] + a[c1])/3.0 - w1);
        const Scalar vol2 = 0.5*Detail::cross(as - bs, a - bs)[c1]*((b[c1] + 2.0*a[c1])/3.0 - w1);
        return abs(vol1 + vol2);
    }

    Scalar volume = 1.0;
    for (const int w : wallIds)
        volume *= spheres[0][Walls<Scalar>::axis(w)] - walls.position(w);
    return abs(volume);
}

//! Volume of the sphere sectors in a tetrahedron, wall vertices contribute no solid
template<class Scalar>
Scalar cellSolidVolume(const std::array<Vertex<Scalar>, 4>& v)
{
    Scalar solid = 0.0;
    for (int i = 0; i < 4; ++i)
        if (v[i].wall < 0)
            solid += solidAngle(v[i].x, v[(i+1)%4].x, v[(i+2)%4].x, v[(i+3)%4].x)*v[i].r*v[i].r*v[i].r/3.0;
    return solid;
}

//! Largest sphere in the void of a tetrahedron with wall vertices
template<class Scalar>
Ball<Scalar> cellInscribedSphere(const std::array<Vertex<Scalar>, 4>& v, const Walls<Scalar>& walls)
{
    Point<Scalar> origin(0.0);
    const std::array<Point<Scalar>, 3> basis{Point<Scalar>{1.0, 0.0, 0.0}, Point<Scalar>{0.0, 1.0, 0.0}, Point<Scalar>{0.0, 0.0, 1.0}};
    return tangentBall<Scalar, 3>(std::vector<Vertex<Scalar>>(v.begin(), v.end()), walls, origin, basis);
}

/*!
 * \brief Geometry of a facet and of the region between it and the power centres p1, p2 of its two
 *        tetrahedra, with up to two wall vertices (YADE FlowEngine)
 *
 * Without walls, this is the facet of facetFluidArea and throatRegion. With one wall, the facet is the
 * trapezoid between the segment of the two sphere centres and its projection onto the wall; with two
 * walls, the rectangle between the sphere centre and its projections onto both walls.
 */
template<class Scalar>
struct FacetGeometry
{
    Point<Scalar> areaVector; //!< oriented to the side of p1
    Scalar fluidArea = 0.0;
    Scalar voidVolume = 0.0;
    std::array<Scalar, 3> solidSurface = {}; //!< zero for walls
    std::array<Scalar, 3> crossSection = {}; //!< circle sector of each sphere in the facet plane, zero for walls
    Scalar hydraulicRadius = 0.0;
    Scalar inscribedRadius = 0.0;
    int numWalls = 0;
};

template<class Scalar>
FacetGeometry<Scalar> facetGeometry(const std::array<Vertex<Scalar>, 3>& v, const Point<Scalar>& p1,
                                    const Point<Scalar>& p2, const Walls<Scalar>& walls)
{
    using std::abs; using std::sqrt;
    FacetGeometry<Scalar> f;
    std::vector<int> real, wallVertices;
    for (int k = 0; k < 3; ++k)
        (v[k].wall < 0 ? real : wallVertices).push_back(k);
    f.numWalls = wallVertices.size();

    const auto sector = [&](int k, const Point<Scalar>& a) {
        const Scalar omega = solidAngle(v[k].x, a, p1, p2);
        return std::array<Scalar, 2>{omega*v[k].r*v[k].r*v[k].r/3.0, omega*v[k].r*v[k].r};
    };

    Scalar solidVolume = 0.0;
    if (f.numWalls == 0)
    {
        f.areaVector = 0.5*Detail::cross(v[0].x - v[2].x, v[1].x - v[2].x);
        for (int i = 0; i < 3; ++i)
            for (int k = 1; k < 3; ++k)
            {
                const auto [volume, surface] = sector(i, v[(i+k)%3].x);
                solidVolume += volume;
                f.solidSurface[i] += surface;
            }
    }
    else if (f.numWalls == 1)
    {
        const int w = v[wallVertices[0]].wall, c = Walls<Scalar>::axis(w);
        const int i2 = real[0], i3 = real[1];
        const auto meanHeight = (walls.position(w) - 0.5*(v[i3].x[c] + v[i2].x[c]))*walls.inwardNormal(w);
        f.areaVector = Detail::cross(meanHeight, v[i3].x - v[i2].x);
        const auto a = walls.project(v[i2].x, w), b = walls.project(v[i3].x, w);
        solidVolume = sector(i2, a)[0] + sector(i2, v[i3].x)[0] + sector(i3, b)[0] + sector(i3, v[i2].x)[0];
        f.solidSurface[i2] = sector(i2, v[wallVertices[0]].x)[1] + sector(i2, v[i3].x)[1];
        f.solidSurface[i3] = sector(i3, v[i2].x)[1] + sector(i3, v[wallVertices[0]].x)[1];
    }
    else
    {
        const int i3 = real[0];
        const auto a = walls.project(v[i3].x, v[wallVertices[0]].wall);
        const auto b = walls.project(v[i3].x, v[wallVertices[1]].wall);
        f.areaVector = Detail::cross(v[i3].x - a, v[i3].x - b);
        solidVolume = sector(i3, a)[0] + sector(i3, b)[0];
        f.solidSurface[i3] = sector(i3, a)[1] + sector(i3, b)[1];
    }
    if (f.areaVector*(p2 - p1) > 0.0)
        f.areaVector *= -1.0;

    f.voidVolume = abs(f.areaVector*(p1 - p2))/3.0 - solidVolume;
    const Scalar solidSurface = f.solidSurface[0] + f.solidSurface[1] + f.solidSurface[2];
    f.hydraulicRadius = f.voidVolume/solidSurface;
    if (f.numWalls == 1)
        f.hydraulicRadius *= 1.0/std::pow(2.0, 0.25);
    else if (f.numWalls == 2)
        f.hydraulicRadius *= 1.0/std::pow(4.0, 0.25);

    // cross-sections with the angles to the actual vertex points, overlap correction as in YADE
    Scalar area = f.areaVector.two_norm();
    for (const int k : real)
    {
        f.crossSection[k] = 0.5*v[k].r*v[k].r*planeAngle(v[k].x, v[(k+1)%3].x, v[(k+2)%3].x);
        area -= f.crossSection[k];
    }
    Scalar overlap = 0.0;
    for (int k = 0; k < 3 && overlap == 0.0; ++k)
        overlap = circularSegmentBeyondEdge(v[k].x, v[k].r, v[(k+1)%3].x, v[(k+2)%3].x);
    f.fluidArea = abs(area + overlap);

    // constriction: largest circle in the facet plane between the sphere sections and the wall lines
    if (f.numWalls == 0)
        f.inscribedRadius = inscribedCircle(std::array<Point<Scalar>, 3>{v[0].x, v[1].x, v[2].x},
                                            std::array<Scalar, 3>{v[0].r, v[1].r, v[2].r}).radius;
    else
    {
        auto e1 = f.areaVector;
        e1 /= e1.two_norm();
        Point<Scalar> u = walls.inwardNormal(v[wallVertices[0]].wall);
        auto e2 = Detail::cross(e1, u);
        e2 /= e2.two_norm();
        const std::array<Point<Scalar>, 2> basis{u, e2};
        std::vector<Vertex<Scalar>> vertices(v.begin(), v.end());
        f.inscribedRadius = tangentBall<Scalar, 2>(vertices, walls, v[real[0]].x, basis).radius;
    }
    return f;
}

} // end namespace Dumux::PoreNetwork::SpherePacking

#endif
