// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup VETest
 * \brief Unit test for column mapping.
 */

#include <config.h>

#include <array>
#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/porousmediumflow/2pve/columnmapping.hh>

int main(int argc, char** argv)
{
    Dumux::initialize(argc, argv);

    static constexpr int dim = 2;
    using Scalar = double;
    using Grid = Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<Scalar, dim>>;
    using GridGeometry = Dumux::CCTpfaFVGridGeometry<typename Grid::LeafGridView>;

    const Dune::FieldVector<Scalar, dim> lowerLeft({1.0, -2.0});
    const Dune::FieldVector<Scalar, dim> upperRight({5.0, 6.0});
    const std::array<int, dim> coarseCells({2, 1});
    const std::array<int, dim> fineCells({2, 4});

    Grid coarseGrid(lowerLeft, upperRight, coarseCells);
    Grid fineGrid(lowerLeft, upperRight, fineCells);

    auto coarseGridGeometry = std::make_shared<GridGeometry>(coarseGrid.leafGridView());
    auto fineGridGeometry = std::make_shared<GridGeometry>(fineGrid.leafGridView());
    const Dumux::VEColumnMapping<GridGeometry, Scalar> mapping(coarseGridGeometry, fineGridGeometry);

    if (mapping.numberOfColumns() != coarseCells[0])
        DUNE_THROW(Dune::Exception, "Expected " << coarseCells[0] << " columns, obtained " << mapping.numberOfColumns());

    for (std::size_t coarseIdx = 0; coarseIdx < mapping.numberOfColumns(); ++coarseIdx)
    {
        const auto& column = mapping.column(coarseIdx);
        if (column.size() != fineCells[1])
            DUNE_THROW(Dune::Exception, "Expected " << fineCells[1] << " fine elements in column " << coarseIdx << ", obtained " << column.size());

        for (std::size_t i = 0; i < column.size(); ++i)
        {
            const auto fineIdx = fineGridGeometry->elementMapper().index(column[i]);
            if (mapping.coarseIndex(fineIdx) != coarseIdx)
                DUNE_THROW(Dune::Exception, "Coarse-to-fine mapping is not invertible for fine element " << fineIdx);

            if (i > 0 && !(column[i-1].geometry().center()[dim-1] < column[i].geometry().center()[dim-1]))
                DUNE_THROW(Dune::Exception, "Column " << coarseIdx << " is not ordered vertically");
        }
    }

    return 0;
}
