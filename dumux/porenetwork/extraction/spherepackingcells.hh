// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkModels
 * \brief Angles of the tetrahedra of a triangulation of a sphere packing at its edges and vertices
 *
 * The dihedral angles of the tetrahedra around an edge between two spheres split the space around
 * the contact of the two spheres among the tetrahedra, and the solid angle of a tetrahedron at a
 * sphere is the part of the sphere surface it contains (R^2 times the solid angle), consistent with
 * the solid volumes of the pores (sphere sectors, see spherepackinggeometry.hh).
 */
#ifndef DUMUX_PNM_EXTRACTION_SPHERE_PACKING_CELLS_HH
#define DUMUX_PNM_EXTRACTION_SPHERE_PACKING_CELLS_HH

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <utility>
#include <vector>

#include <dumux/porenetwork/extraction/spherepackinggeometry.hh>
#include <dumux/porenetwork/extraction/spherepackingnetwork.hh>
#include <dumux/porenetwork/extraction/spherepackingwalls.hh>

namespace Dumux::PoreNetwork::SpherePacking {

//! Positions of the vertices of a triangulation with walls: the sphere centres, then the wall vertices
template<class Scalar>
std::vector<Point<Scalar>> vertexPositions(const std::vector<Point<Scalar>>& centers, const Walls<Scalar>& walls)
{
    auto positions = centers;
    for (int k = 0; k < 6; ++k)
        positions.push_back(walls.center(k));
    return positions;
}

//! Interior angle of the tetrahedron (a, b, c, d) at its edge (a, b)
template<class Scalar>
Scalar dihedralAngle(const Point<Scalar>& a, const Point<Scalar>& b, const Point<Scalar>& c, const Point<Scalar>& d)
{
    using std::atan2;
    auto e = b - a;
    e /= e.two_norm();
    auto u = c - a, v = d - a;
    u.axpy(-(u*e), e);
    v.axpy(-(v*e), e);
    return atan2(Detail::cross(u, v).two_norm(), u*v);
}

/*!
 * \brief Edge of a triangulation between two spheres, with the tetrahedra around it and their
 *        dihedral angles at the edge
 */
template<class Scalar>
struct SphereEdge
{
    std::array<int, 2> spheres; //!< sphere indices, spheres[0] < spheres[1]
    std::vector<std::size_t> cells;
    std::vector<Scalar> dihedralAngles;
};

/*!
 * \brief All edges between two spheres of a triangulation whose vertices numSpheres and above are
 *        walls (or other non-sphere vertices), sorted by sphere indices
 *
 * The dihedral angles around an edge in the interior of the triangulated region add up to 2 pi.
 */
template<class Scalar>
std::vector<SphereEdge<Scalar>> sphereEdges(const Triangulation& triangulation,
                                            const std::vector<Point<Scalar>>& positions,
                                            int numSpheres)
{
    std::map<std::pair<int, int>, std::size_t> index;
    std::vector<SphereEdge<Scalar>> edges;
    for (std::size_t c = 0; c < triangulation.cells.size(); ++c)
    {
        const auto& v = triangulation.cells[c];
        for (int a = 0; a < 4; ++a)
            for (int b = a + 1; b < 4; ++b)
            {
                if (v[a] >= numSpheres || v[b] >= numSpheres)
                    continue;
                int others[2], m = 0;
                for (int k = 0; k < 4; ++k)
                    if (k != a && k != b)
                        others[m++] = v[k];
                const auto key = std::minmax(v[a], v[b]);
                auto [it, inserted] = index.try_emplace(key, edges.size());
                if (inserted)
                    edges.push_back({{key.first, key.second}, {}, {}});
                auto& edge = edges[it->second];
                edge.cells.push_back(c);
                edge.dihedralAngles.push_back(dihedralAngle(positions[v[a]], positions[v[b]],
                                                            positions[others[0]], positions[others[1]]));
            }
    }
    std::sort(edges.begin(), edges.end(), [](const auto& e0, const auto& e1) { return e0.spheres < e1.spheres; });
    return edges;
}

/*!
 * \brief Solid angle of every tetrahedron at each of its four vertices
 *
 * For a sphere vertex, R^2 times this angle is the part of the sphere surface inside the tetrahedron.
 */
template<class Scalar>
std::vector<std::array<Scalar, 4>> cellSolidAngles(const Triangulation& triangulation,
                                                   const std::vector<Point<Scalar>>& positions)
{
    std::vector<std::array<Scalar, 4>> angles(triangulation.cells.size());
    for (std::size_t c = 0; c < triangulation.cells.size(); ++c)
    {
        const auto& v = triangulation.cells[c];
        for (int i = 0; i < 4; ++i)
            angles[c][i] = solidAngle(positions[v[i]], positions[v[(i+1)%4]], positions[v[(i+2)%4]], positions[v[(i+3)%4]]);
    }
    return angles;
}

} // end namespace Dumux::PoreNetwork::SpherePacking

#endif
