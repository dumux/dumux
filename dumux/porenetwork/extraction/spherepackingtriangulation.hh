// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkModels
 * \brief Regular triangulation of a sphere packing (requires CGAL)
 */
#ifndef DUMUX_PNM_EXTRACTION_SPHERE_PACKING_TRIANGULATION_HH
#define DUMUX_PNM_EXTRACTION_SPHERE_PACKING_TRIANGULATION_HH

#include <cstddef>
#include <utility>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Regular_triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include <dumux/porenetwork/extraction/spherepackingnetwork.hh>

namespace Dumux::PoreNetwork::SpherePacking {

namespace Detail {

//! Regular triangulation of weighted points (x, w = r^2) with their indices as vertex info
template<class Point, class Scalar>
Triangulation regularTriangulation(const std::vector<Point>& centers, const std::vector<Scalar>& radii,
                                   std::size_t numSpheres)
{
    using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
    using VertexBase = CGAL::Triangulation_vertex_base_with_info_3<int, Kernel, CGAL::Regular_triangulation_vertex_base_3<Kernel>>;
    using CellBase = CGAL::Triangulation_cell_base_with_info_3<int, Kernel, CGAL::Regular_triangulation_cell_base_3<Kernel>>;
    using DataStructure = CGAL::Triangulation_data_structure_3<VertexBase, CellBase>;
    using RegularTriangulation = CGAL::Regular_triangulation_3<Kernel, DataStructure>;
    using WeightedPoint = Kernel::Weighted_point_3;

    std::vector<std::pair<WeightedPoint, int>> points;
    points.reserve(centers.size());
    for (std::size_t i = 0; i < centers.size(); ++i)
    {
        const auto& x = centers[i];
        points.emplace_back(WeightedPoint(Kernel::Point_3(x[0], x[1], x[2]), radii[i]*radii[i]), static_cast<int>(i));
    }

    RegularTriangulation rt(points.begin(), points.end());

    Triangulation result;
    std::size_t numSphereVertices = 0;
    for (auto v = rt.finite_vertices_begin(); v != rt.finite_vertices_end(); ++v)
        numSphereVertices += static_cast<std::size_t>(v->info()) < numSpheres;
    result.numHiddenSpheres = numSpheres - numSphereVertices;

    int index = 0;
    for (auto cell = rt.finite_cells_begin(); cell != rt.finite_cells_end(); ++cell)
        cell->info() = index++;

    result.cells.resize(index);
    result.neighbors.resize(index);
    for (auto cell = rt.finite_cells_begin(); cell != rt.finite_cells_end(); ++cell)
    {
        const int c = cell->info();
        for (int j = 0; j < 4; ++j)
        {
            result.cells[c][j] = cell->vertex(j)->info();
            const auto neighbor = cell->neighbor(j);
            result.neighbors[c][j] = rt.is_infinite(neighbor) ? -1 : neighbor->info();
        }
    }
    return result;
}

} // end namespace Detail

/*!
 * \brief Regular (weighted Delaunay) triangulation of the sphere centres with weights r^2
 *
 * Its dual is the power diagram (radical Voronoi tessellation) of the spheres.
 * Exact predicates make the combinatorics robust; degenerate configurations such as lattices
 * are resolved by symbolic perturbation.
 */
template<class Point, class Scalar>
Triangulation regularTriangulation(const std::vector<Point>& centers, const std::vector<Scalar>& radii)
{ return Detail::regularTriangulation(centers, radii, centers.size()); }

/*!
 * \brief Regular triangulation of the spheres and the six walls of a box
 *
 * The walls are the vertices numSpheres + k, spheres of very large radius whose surface is the wall
 * plane (see Walls).
 */
template<class Point, class Scalar>
Triangulation regularTriangulation(const std::vector<Point>& centers, const std::vector<Scalar>& radii,
                                   const Walls<Scalar>& walls)
{
    auto points = centers;
    auto weights = radii;
    for (int k = 0; k < 6; ++k)
    {
        points.push_back(walls.center(k));
        weights.push_back(walls.radius());
    }
    return Detail::regularTriangulation(points, weights, centers.size());
}

} // end namespace Dumux::PoreNetwork::SpherePacking

#endif
