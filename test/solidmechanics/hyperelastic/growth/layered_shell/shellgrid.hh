// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Structured cube grids of a quarter of a two-layer annulus (2D) and of an eighth of a
 *        two-layer spherical shell (3D)
 *
 * The radial nodes are uniform in each layer, so the layer interface is a grid surface. In 3D the
 * octant is covered by the three cube faces \f$ x = 1 \f$, \f$ y = 1 \f$, \f$ z = 1 \f$ of the
 * equiangular cubed sphere, each with \f$ m \times m \f$ cells.
 */
#ifndef DUMUX_TEST_HYPERELASTIC_GROWTH_LAYERED_SHELL_GRID_HH
#define DUMUX_TEST_HYPERELASTIC_GROWTH_LAYERED_SHELL_GRID_HH

#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/common/gridfactory.hh>

namespace Dumux {

template<class Scalar>
std::vector<Scalar> layeredShellRadialNodes(const std::array<Scalar, 3>& radii, int cellsInner, int cellsOuter)
{
    std::vector<Scalar> nodes;
    for (int i = 0; i < cellsInner; ++i)
        nodes.push_back(radii[0] + (radii[1] - radii[0])*i/cellsInner);
    for (int i = 0; i <= cellsOuter; ++i)
        nodes.push_back(radii[1] + (radii[2] - radii[1])*i/cellsOuter);
    return nodes;
}

/*!
 * \brief Creates the grid
 * \param radii inner radius, interface radius and outer radius
 * \param cells radial cells of the inner layer, radial cells of the outer layer and angular cells
 *        (2D: over the quarter circle, 3D: along each edge of a cubed-sphere face)
 */
template<class Grid, class Scalar>
std::unique_ptr<Grid> createLayeredShellGrid(const std::array<Scalar, 3>& radii, const std::array<int, 3>& cells)
{
    static constexpr int dim = Grid::dimension;
    using Position = Dune::FieldVector<typename Grid::ctype, dim>;
    const auto rNodes = layeredShellRadialNodes(radii, cells[0], cells[1]);
    const int nr = rNodes.size();
    const int m = cells[2];

    Dune::GridFactory<Grid> factory;
    std::map<std::array<long long, dim>, unsigned int> vertexIndex;
    const auto insertVertex = [&](const Position& pos)
    {
        std::array<long long, dim> key;
        for (int i = 0; i < dim; ++i)
            key[i] = std::llround(pos[i]*1e10);
        const auto [it, inserted] = vertexIndex.try_emplace(key, vertexIndex.size());
        if (inserted)
            factory.insertVertex(pos);
        return it->second;
    };

    // inserts the cell with corners corner(a, b, c), a, b, c in {0, 1}, positively oriented
    const auto insertCell = [&](const auto& corner)
    {
        std::vector<Position> pos;
        std::vector<unsigned int> indices;
        for (int c = 0; c < (dim == 3 ? 2 : 1); ++c)
            for (int b = 0; b < 2; ++b)
                for (int a = 0; a < 2; ++a)
                {
                    pos.push_back(corner(a, b, c));
                    indices.push_back(insertVertex(pos.back()));
                }

        // orientation from the edges at the first corner
        Dune::FieldMatrix<Scalar, dim, dim> J;
        for (int i = 0; i < dim; ++i)
            J[i] = pos[1 << i] - pos[0];
        if (J.determinant() < 0.0)
            for (int k = 0; k < static_cast<int>(indices.size()); k += 2)
                std::swap(indices[k], indices[k+1]);

        factory.insertElement(Dune::GeometryTypes::cube(dim), indices);
    };

    if constexpr (dim == 2)
    {
        const auto point = [&](int i, int j) {
            const Scalar phi = 0.5*M_PI*j/m;
            return Position{rNodes[i]*std::cos(phi), rNodes[i]*std::sin(phi)};
        };
        for (int i = 0; i < nr - 1; ++i)
            for (int j = 0; j < m; ++j)
                insertCell([&](int a, int b, int) { return point(i + a, j + b); });
    }
    else
    {
        // point on face f of the cubed sphere (the coordinate f is 1) at angular indices (j, k)
        const auto point = [&](int f, int i, int j, int k) {
            Position d;
            d[f] = 1.0;
            d[(f + 1) % 3] = std::tan(0.25*M_PI*j/m);
            d[(f + 2) % 3] = std::tan(0.25*M_PI*k/m);
            d *= rNodes[i]/d.two_norm();
            return d;
        };
        for (int f = 0; f < 3; ++f)
            for (int i = 0; i < nr - 1; ++i)
                for (int j = 0; j < m; ++j)
                    for (int k = 0; k < m; ++k)
                        insertCell([&](int a, int b, int c) { return point(f, i + c, j + a, k + b); });
    }

    return std::unique_ptr<Grid>(factory.createGrid());
}

} // end namespace Dumux

#endif
