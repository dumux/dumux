// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test that the finite element Jacobian pattern is unchanged for bases over a grid.
 */
#include <config.h>

#include <vector>
#include <memory>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/istl/matrixindexset.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/concepts/functionspacebasis_.hh>
#include <dumux/assembly/jacobianpattern.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>
#include <dumux/discretization/projection/l2_projection.hh>

int main(int argc, char** argv)
{
    using namespace Dumux;
    initialize(argc, argv);

    using Grid = Dune::YaspGrid<2>;
    Grid grid{{1.0, 1.0}, {4, 4}};
    const auto gridView = grid.leafGridView();

    using GridGeometry = BoxFVGridGeometry<double, Grid::LeafGridView>;
    const auto gridGeometry = std::make_shared<GridGeometry>(gridView);

    const FEBasisFromCVFEGridDiscretization<GridGeometry> basis{*gridGeometry};
    static_assert(Concept::ProjectionBasis<FEBasisFromCVFEGridDiscretization<GridGeometry>>);
    static_assert(Concept::EntityRangeProvider<FEBasisFromCVFEGridDiscretization<GridGeometry>>);

    if (basis.size() != 25)
        DUNE_THROW(Dune::InvalidStateException,
                   "Expected 25 vertices on a 4x4 grid, got " << basis.size());

    const auto pattern = getFEJacobianPattern(basis);

    // independently: two vertices couple exactly when they share an element
    Dune::MatrixIndexSet expected;
    expected.resize(basis.size(), basis.size());
    for (const auto& element : elements(gridView))
    {
        const auto numVertices = element.subEntities(Grid::dimension);
        for (unsigned int i = 0; i < numVertices; ++i)
            for (unsigned int j = 0; j < numVertices; ++j)
                expected.add(gridGeometry->vertexMapper().subIndex(element, i, Grid::dimension),
                             gridGeometry->vertexMapper().subIndex(element, j, Grid::dimension));
    }

    if (pattern.rows() != expected.rows())
        DUNE_THROW(Dune::InvalidStateException, "Pattern row count mismatch");

    for (std::size_t row = 0; row < pattern.rows(); ++row)
        if (pattern.rowsize(row) != expected.rowsize(row))
            DUNE_THROW(Dune::InvalidStateException,
                       "Row " << row << " has " << pattern.rowsize(row)
                              << " entries, expected " << expected.rowsize(row));

    std::cout << "FE Jacobian pattern over a grid view: " << pattern.rows()
              << " rows, " << pattern.size() << " entries, matches the element-wise construction"
              << std::endl;
    return 0;
}
