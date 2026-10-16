// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Linear
 * \brief Prolongation from the P1/Q1 vertex space to a higher-order Lagrange space
 */
#ifndef DUMUX_LINEAR_LAGRANGE_PROLONGATION_HH
#define DUMUX_LINEAR_LAGRANGE_PROLONGATION_HH

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/localfunctions/lagrange/lagrangelfecache.hh>

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief Assemble the prolongation matrix from the P1/Q1 space on the grid vertices
 *        to the Lagrange space of a grid geometry.
 *
 * The coarse space is the space of piecewise linear (simplices) or multilinear (cubes)
 * functions on the same grid, with one coarse dof per grid vertex (numbered by the
 * vertex mapper of the grid view). Row \f$i\f$ of the prolongation \f$P\f$ contains the
 * values of the P1/Q1 basis functions of an element at the position of fine dof \f$i\f$,
 * i.e. \f$P\f$ interpolates a P1/Q1 function in the nodes of the fine space.
 * Works for all grid geometries with the local dof interface (FE and CVFE hybrid PQk
 * grid geometries). The grid geometry is only used to iterate over the local dofs
 * of the elements, the vertex dofs of the fine space need not be identified.
 *
 * \tparam blockSize number of equations per dof. The entries of the returned matrix are
 *                   blockSize x blockSize multiples of the identity such that the
 *                   result is applicable to systems with blockSize unknowns per dof.
 * \tparam Scalar the scalar type
 * \param gridGeometry the grid geometry of the fine space
 * \param tolerance weights below this value are not stored
 * \note Dofs shared between elements get the same weights from all elements, which is
 *       checked up to a tolerance of 1e-10.
 */
template<int blockSize = 1, class Scalar = double, class GridGeometry>
Dune::BCRSMatrix<Dune::FieldMatrix<Scalar, blockSize, blockSize>>
makeLagrangeP1Prolongation(const GridGeometry& gridGeometry, Scalar tolerance = 1e-12)
{
    using GridView = typename GridGeometry::GridView;
    using Block = Dune::FieldMatrix<Scalar, blockSize, blockSize>;
    using Matrix = Dune::BCRSMatrix<Block>;
    using Row = std::vector<std::pair<std::size_t, Scalar>>;
    constexpr int dim = GridView::dimension;

    const auto& gridView = gridGeometry.gridView();
    Dune::MultipleCodimMultipleGeomTypeMapper<GridView> vertexMapper(gridView, Dune::mcmgVertexLayout());
    Dune::LagrangeLocalFiniteElementCache<typename GridView::ctype, Scalar, dim, 1> p1Cache;

    std::vector<Row> rows(gridGeometry.numDofs());
    std::vector<bool> isSet(gridGeometry.numDofs(), false);
    std::vector<Dune::FieldVector<Scalar, 1>> shapeValues;

    auto fvGeometry = localView(gridGeometry);
    for (const auto& element : elements(gridView))
    {
        fvGeometry.bind(element);
        const auto& p1Basis = p1Cache.get(element.type()).localBasis();
        for (const auto& localDof : localDofs(fvGeometry))
        {
            p1Basis.evaluateFunction(ipData(fvGeometry, localDof).local(), shapeValues);

            Row row;
            for (std::size_t i = 0; i < shapeValues.size(); ++i)
                if (std::abs(shapeValues[i][0]) > tolerance)
                    row.emplace_back(vertexMapper.subIndex(element, i, dim), shapeValues[i][0]);
            std::ranges::sort(row);

            const auto dofIdx = localDof.dofIndex();
            if (!isSet[dofIdx])
            {
                rows[dofIdx] = std::move(row);
                isSet[dofIdx] = true;
            }
            else
            {
                const auto& ref = rows[dofIdx];
                bool consistent = (ref.size() == row.size());
                for (std::size_t k = 0; consistent && k < ref.size(); ++k)
                    consistent = ref[k].first == row[k].first && std::abs(ref[k].second - row[k].second) < 1e-10;
                if (!consistent)
                    DUNE_THROW(Dune::InvalidStateException,
                               "Inconsistent prolongation weights for dof " << dofIdx);
            }
        }
    }

    Matrix P(gridGeometry.numDofs(), vertexMapper.size(), Matrix::random);
    for (std::size_t i = 0; i < rows.size(); ++i)
        P.setrowsize(i, rows[i].size());
    P.endrowsizes();
    for (std::size_t i = 0; i < rows.size(); ++i)
        for (const auto& [col, w] : rows[i])
            P.addindex(i, col);
    P.endindices();
    for (std::size_t i = 0; i < rows.size(); ++i)
        for (const auto& [col, w] : rows[i])
        {
            Block b(0.0);
            for (int k = 0; k < blockSize; ++k)
                b[k][k] = w;
            P[i][col] = b;
        }
    return P;
}

} // end namespace Dumux

#endif
