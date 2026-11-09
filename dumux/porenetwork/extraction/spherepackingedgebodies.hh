// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkModels
 * \brief Pore network of a sphere packing with axisymmetric solid bodies at the edges between spheres,
 *        e.g. the lenses of overlapping spheres or cement rings at contacts, and with solid tetrahedra
 *
 * A body on the edge between spheres i and j is rotationally symmetric about the line through their centres
 * and lies near that line. Its part in a tetrahedron (or throat region) around the edge is its volume times the
 * dihedral angle of the tetrahedron at the edge over 2 pi, and every facet containing the edge cuts it in a
 * half-plane bounded by the axis. This holds while the body is closer to the edge than the faces of the
 * tetrahedra not containing the edge.
 */
#ifndef DUMUX_PNM_EXTRACTION_SPHERE_PACKING_EDGE_BODIES_HH
#define DUMUX_PNM_EXTRACTION_SPHERE_PACKING_EDGE_BODIES_HH

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <numbers>
#include <queue>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dumux/porenetwork/extraction/spherepackingcells.hh>
#include <dumux/porenetwork/extraction/spherepackinggeometry.hh>
#include <dumux/porenetwork/extraction/spherepackingnetwork.hh>

namespace Dumux::PoreNetwork::SpherePacking {

/*!
 * \brief Axisymmetric solid body at an edge between two spheres
 *
 * The quantities are changes relative to the spheres alone: solid volume added to the void (negative for
 * solid that the sphere sectors count twice), solid area in a half-plane bounded by the axis, and change of
 * the solid surface in contact with the fluid.
 */
template<class Scalar>
struct EdgeBody
{
    Scalar solidVolume = 0.0;
    Scalar halfSection = 0.0;
    Scalar wettedSurfaceChange = 0.0;

    EdgeBody& operator+=(const EdgeBody& other)
    {
        solidVolume += other.solidVolume;
        halfSection += other.halfSection;
        wettedSurfaceChange += other.wettedSurfaceChange;
        return *this;
    }
};

//! Bodies by sphere pair (i, j) with i < j
template<class Scalar>
using EdgeBodies = std::map<std::pair<int, int>, EdgeBody<Scalar>>;

/*!
 * \brief Lens in which two spheres of radii r1, r2 at the given centre distance overlap
 *
 * capHeights are the depths of the two spheres beyond the radical plane, radius the radius of the circle in
 * which the sphere surfaces meet, area the area of the lens in a plane through both centres.
 */
template<class Scalar>
struct Lens
{
    Scalar volume = 0.0;
    Scalar area = 0.0;
    Scalar radius = 0.0;
    Scalar radicalDistance = 0.0; //!< distance of the radical plane from the centre of sphere 1
    std::array<Scalar, 2> capHeights = {};
};

//! Lens of two spheres; zero unless they overlap without one containing the other
template<class Scalar>
Lens<Scalar> sphereLens(Scalar r1, Scalar r2, Scalar distance)
{
    using std::sqrt; using std::acos; using std::abs;
    Lens<Scalar> lens;
    lens.radicalDistance = (distance*distance + r1*r1 - r2*r2)/(2.0*distance);
    if (distance >= r1 + r2 || distance <= abs(r1 - r2))
        return lens;

    const Scalar x1 = lens.radicalDistance, x2 = distance - x1;
    lens.radius = sqrt(std::max(r1*r1 - x1*x1, Scalar(0.0)));
    lens.capHeights = {r1 - x1, r2 - x2};
    const auto& h = lens.capHeights;
    lens.volume = std::numbers::pi/3.0*(h[0]*h[0]*(3.0*r1 - h[0]) + h[1]*h[1]*(3.0*r2 - h[1]));
    lens.area = r1*r1*acos(std::clamp(x1/r1, Scalar(-1.0), Scalar(1.0))) - x1*lens.radius
                + r2*r2*acos(std::clamp(x2/r2, Scalar(-1.0), Scalar(1.0))) - x2*lens.radius;
    return lens;
}

/*!
 * \brief Edge body that removes the lens of two overlapping spheres from the solid counted twice by the sphere
 *        sectors: the lens volume, half its area and the two sphere caps inside the other sphere
 */
template<class Scalar>
EdgeBody<Scalar> overlapCorrection(Scalar r1, Scalar r2, Scalar distance)
{
    const auto lens = sphereLens(r1, r2, distance);
    const Scalar caps = 2.0*std::numbers::pi*(r1*lens.capHeights[0] + r2*lens.capHeights[1]);
    return {-lens.volume, -0.5*lens.area, -caps};
}

/*!
 * \brief Point in the plane of three spheres with equal power with respect to all of them
 */
template<class Scalar>
Point<Scalar> facetPowerCenter(const std::array<Point<Scalar>, 3>& x, const std::array<Scalar, 3>& r)
{
    const auto e1 = x[1] - x[0], e2 = x[2] - x[0];
    // f = x0 + s e1 + t e2 with 2 (x_k - x0).(f - x0) = |x_k - x0|^2 - r_k^2 + r_0^2
    const Scalar a11 = e1*e1, a12 = e1*e2, a22 = e2*e2;
    const Scalar b1 = 0.5*(a11 - r[1]*r[1] + r[0]*r[0]);
    const Scalar b2 = 0.5*(a22 - r[2]*r[2] + r[0]*r[0]);
    const Scalar det = a11*a22 - a12*a12;
    const Scalar s = (b1*a22 - b2*a12)/det;
    const Scalar t = (a11*b2 - a12*b1)/det;
    auto f = x[0];
    f.axpy(s, e1);
    f.axpy(t, e2);
    return f;
}

namespace Detail {

//! Fluid area of a facet as facetFluidArea, but signed, so that sectors counted twice may make it negative
template<class Scalar>
Scalar signedFacetFluidArea(const std::array<Point<Scalar>, 3>& x, const std::array<Scalar, 3>& r)
{
    Scalar area = triangleArea(x[0], x[1], x[2]);
    for (int i = 0; i < 3; ++i)
        area -= 0.5*r[i]*r[i]*planeAngle(x[i], x[(i+1)%3], x[(i+2)%3]);

    for (int i = 0; i < 3; ++i)
        if (const Scalar segment = circularSegmentBeyondEdge(x[i], r[i], x[(i+1)%3], x[(i+2)%3]); segment > 0.0)
            return area + segment;

    return area;
}

} // end namespace Detail

/*!
 * \brief Overlap corrections of all edges of a triangulation between overlapping spheres
 *
 * Vertices numSpheres and above are walls and get no correction. Throws Dune::InvalidStateException if three
 * spheres of a tetrahedron share a common point, where the lenses overlap and the correction is not exact.
 */
template<class Scalar>
EdgeBodies<Scalar> overlapCorrections(const Triangulation& triangulation,
                                      const std::vector<Point<Scalar>>& centers,
                                      const std::vector<Scalar>& radii)
{
    const int numSpheres = static_cast<int>(centers.size());
    EdgeBodies<Scalar> bodies;
    for (const auto& v : triangulation.cells)
    {
        for (int a = 0; a < 4; ++a)
            for (int b = a + 1; b < 4; ++b)
            {
                if (v[a] >= numSpheres || v[b] >= numSpheres)
                    continue;
                const auto key = std::minmax(v[a], v[b]);
                const Scalar distance = (centers[key.second] - centers[key.first]).two_norm();
                if (distance < radii[key.first] + radii[key.second] && !bodies.count(key))
                    bodies[key] = overlapCorrection(radii[key.first], radii[key.second], distance);
            }

        // a point inside three spheres is a facet point of negative power with respect to all three
        for (int skip = 0; skip < 4; ++skip)
        {
            std::array<Point<Scalar>, 3> x;
            std::array<Scalar, 3> r;
            bool spheresOnly = true;
            for (int k = 0, m = 0; k < 4; ++k)
            {
                if (k == skip)
                    continue;
                spheresOnly = spheresOnly && v[k] < numSpheres;
                if (spheresOnly)
                {
                    x[m] = centers[v[k]];
                    r[m] = radii[v[k]];
                }
                ++m;
            }
            if (!spheresOnly)
                continue;
            const auto f = facetPowerCenter(x, r);
            const Scalar power = (f - x[0]).two_norm2() - r[0]*r[0];
            if (power < 0.0)
                DUNE_THROW(Dune::InvalidStateException, "Three spheres overlap in a common region");
        }
    }
    return bodies;
}

/*!
 * \brief Pore network of the spheres in a region (see extractNetwork with an ExtractionRegion), with edge
 *        bodies and with solid tetrahedra
 *
 * Every pore loses the parts of the bodies at its edges, every facet their half-sections, every throat region
 * their parts in it from its void volume and their change of the wetted surface from its solid surface.
 * solidCells[c] marks tetrahedron c as solid: it is no pore and its facets are closed throats; empty means none.
 * Tetrahedra with a vertex numSpheres = centers.size() or above (walls) are left out, so the triangulation may
 * include walls.
 */
template<class Scalar>
Network<Scalar> extractNetwork(const Triangulation& triangulation,
                               const std::vector<Point<Scalar>>& centers,
                               const std::vector<Scalar>& radii,
                               const ExtractionRegion<Scalar>& region,
                               const EdgeBodies<Scalar>& bodies,
                               const std::vector<char>& solidCells = {})
{
    using std::abs;
    constexpr Scalar twoPi = 2.0*std::numbers::pi;
    const int numSpheres = static_cast<int>(centers.size());
    const std::size_t numCells = triangulation.cells.size();
    const auto body = [&](int i, int j) -> const EdgeBody<Scalar>* {
        const auto it = bodies.find(std::minmax(i, j));
        return it == bodies.end() ? nullptr : &it->second;
    };
    const auto isSolid = [&](std::size_t c) { return !solidCells.empty() && solidCells[c]; };

    std::vector<Detail::CellRecord<Scalar>> cells(numCells);
    std::size_t numSolidFacets = 0;
    for (std::size_t c = 0; c < numCells; ++c)
    {
        const auto& v = triangulation.cells[c];
        if (std::any_of(v.begin(), v.end(), [&](int vertex) { return vertex >= numSpheres; }) || isSolid(c))
            continue;
        std::array<Point<Scalar>, 4> x;
        std::array<Scalar, 4> r;
        for (int j = 0; j < 4; ++j)
        {
            x[j] = centers[v[j]];
            r[j] = radii[v[j]];
        }
        const auto centroid = (x[0] + x[1] + x[2] + x[3])/4.0;
        bool inside = true;
        for (int axis = 0; axis < 3; ++axis)
            inside = inside && centroid[axis] >= region.lower[axis] && centroid[axis] <= region.upper[axis];
        if (!inside)
            continue;

        auto& cell = cells[c];
        cell.kept = true;
        cell.pore.position = powerCenter(x, r);
        cell.pore.bulkVolume = tetrahedronVolume(x);
        cell.pore.volume = cell.pore.bulkVolume - tetrahedronSolidVolume(x, r);
        for (int a = 0; a < 4; ++a)
            for (int b = a + 1; b < 4; ++b)
                if (const auto* e = body(v[a], v[b]))
                {
                    int others[2], m = 0;
                    for (int k = 0; k < 4; ++k)
                        if (k != a && k != b)
                            others[m++] = k;
                    cell.pore.volume -= e->solidVolume*dihedralAngle(x[a], x[b], x[others[0]], x[others[1]])/twoPi;
                }
        cell.pore.inscribedRadius = inscribedSphere(x, r).radius;
        if (cell.pore.inscribedRadius <= 0.0)
            cell.pore.inscribedRadius = clearance(cell.pore.position, x, r);
        for (int face = 0; face < 6; ++face)
        {
            const int axis = face/2;
            for (const auto& p : x)
            {
                const Scalar depth = face % 2 == 0 ? region.lower[axis] - p[axis] : p[axis] - region.upper[axis];
                if (depth > cell.labelPriority)
                {
                    cell.labelPriority = depth;
                    cell.pore.label = region.boundaryLabels[face];
                }
            }
        }
    }

    std::vector<Detail::FacetRecord<Scalar>> facets;
    for (std::size_t c = 0; c < numCells; ++c)
    {
        if (!cells[c].kept)
            continue;
        for (int j = 0; j < 4; ++j)
        {
            const int n = triangulation.neighbors[c][j];
            if (n >= 0 && isSolid(n))
                ++numSolidFacets;
            if (n < 0 || !cells[n].kept || static_cast<std::size_t>(n) < c)
                continue;

            Detail::FacetRecord<Scalar> facet;
            facet.cells = {c, static_cast<std::size_t>(n)};
            std::array<Point<Scalar>, 3> x;
            std::array<Scalar, 3> r;
            for (int k = 0, m = 0; k < 4; ++k)
            {
                if (k == j)
                    continue;
                facet.vertices[m] = triangulation.cells[c][k];
                x[m] = centers[facet.vertices[m]];
                r[m] = radii[facet.vertices[m]];
                ++m;
            }

            const auto& p1 = cells[c].pore.position;
            const auto& p2 = cells[n].pore.position;
            facet.length = (p1 - p2).two_norm();
            facet.fluidArea = Detail::signedFacetFluidArea(x, r);
            facet.inscribedRadius = inscribedCircle(x, r).radius;
            auto throat = throatRegion(x, r, p1, p2);
            for (int k = 0; k < 3; ++k)
                if (const auto* e = body(facet.vertices[k], facet.vertices[(k+1)%3]))
                {
                    const Scalar share = dihedralAngle(x[k], x[(k+1)%3], p1, p2)/twoPi;
                    facet.fluidArea -= e->halfSection;
                    throat.voidVolume -= e->solidVolume*share;
                    throat.solidSurface += e->wettedSurfaceChange*share;
                }
            facet.hydraulicRadius = facet.fluidArea > 0.0 ? throat.hydraulicRadius() : 0.0;
            Scalar wettedPerimeter = 0.0;
            for (int k = 0; k < 3; ++k)
                wettedPerimeter += r[k]*planeAngle(x[k], x[(k+1)%3], x[(k+2)%3]);
            facet.shapeFactor = facet.fluidArea/(wettedPerimeter*wettedPerimeter);
            const auto force = facetForce(x, r, p1, p2);
            facet.normal = force.normal;
            facet.forceWeights = force.weights;
            facets.push_back(facet);
        }
    }

    auto network = Detail::assembleNetwork(cells, facets, std::vector<Detail::WallRecord<Scalar>>{},
                                           region.mergeDistance, false, region.boundaryLabels);
    network.numClosedThroats += numSolidFacets;
    return network;
}

/*!
 * \brief The network without the pores that no chain of throats connects to a pore with one of the labels
 *
 * Throats, facets, wall pressures and the pores of the tetrahedra are renumbered accordingly.
 */
template<class Scalar>
Network<Scalar> connectedNetwork(const Network<Scalar>& network, const std::vector<int>& labels)
{
    const std::size_t n = network.pores.size();
    std::vector<std::vector<std::size_t>> adjacent(n);
    for (const auto& t : network.throats)
    {
        adjacent[t.pores[0]].push_back(t.pores[1]);
        adjacent[t.pores[1]].push_back(t.pores[0]);
    }
    std::vector<char> reached(n, 0);
    std::queue<std::size_t> front;
    for (std::size_t p = 0; p < n; ++p)
        if (std::find(labels.begin(), labels.end(), network.pores[p].label) != labels.end())
        {
            reached[p] = 1;
            front.push(p);
        }
    while (!front.empty())
    {
        const auto p = front.front();
        front.pop();
        for (const auto q : adjacent[p])
            if (!reached[q])
            {
                reached[q] = 1;
                front.push(q);
            }
    }

    Network<Scalar> result;
    result.numClosedThroats = network.numClosedThroats;
    result.numMergedTetrahedra = network.numMergedTetrahedra;
    result.numIsolatedPores = network.numIsolatedPores + std::count(reached.begin(), reached.end(), 0);
    std::vector<long> index(n, -1);
    for (std::size_t p = 0; p < n; ++p)
        if (reached[p])
        {
            index[p] = result.pores.size();
            result.pores.push_back(network.pores[p]);
        }
    for (auto t : network.throats)
        if (reached[t.pores[0]])
        {
            t.pores = {static_cast<std::size_t>(index[t.pores[0]]), static_cast<std::size_t>(index[t.pores[1]])};
            result.throats.push_back(t);
        }
    for (auto f : network.facets)
        if (reached[f.pores[0]])
        {
            f.pores = {static_cast<std::size_t>(index[f.pores[0]]), static_cast<std::size_t>(index[f.pores[1]])};
            result.facets.push_back(f);
        }
    for (auto w : network.wallPressures)
        if (reached[w.pore])
        {
            w.pore = index[w.pore];
            result.wallPressures.push_back(w);
        }
    result.cellPores.resize(network.cellPores.size());
    for (std::size_t c = 0; c < network.cellPores.size(); ++c)
        result.cellPores[c] = network.cellPores[c] >= 0 ? index[network.cellPores[c]] : -1;
    return result;
}

} // end namespace Dumux::PoreNetwork::SpherePacking

#endif
