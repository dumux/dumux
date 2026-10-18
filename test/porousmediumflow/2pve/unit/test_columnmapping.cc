// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPVETests
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

template<int dim>
void checkMapping()
{
    using Scalar = double;
    using Grid = Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<Scalar, dim>>;
    using GridGeometry = Dumux::CCTpfaFVGridGeometry<typename Grid::LeafGridView>;
    using CoarseGrid = Dune::YaspGrid<dim-1, Dune::EquidistantOffsetCoordinates<Scalar, dim-1>>;
    using CoarseGridGeometry = Dumux::CCTpfaFVGridGeometry<typename CoarseGrid::LeafGridView>;

    Dune::FieldVector<Scalar, dim> lowerLeft(1.0), upperRight(5.0);
    lowerLeft[dim-1] = -2.0;
    upperRight[dim-1] = 6.0;
    Dune::FieldVector<Scalar, dim-1> lowerHorizontal(1.0), upperHorizontal(5.0);
    std::array<int, dim-1> coarseCells;
    coarseCells.fill(2);
    std::array<int, dim> fineCells;
    fineCells.fill(2);
    fineCells[dim-1] = 4;
    CoarseGrid coarseGrid(lowerHorizontal, upperHorizontal, coarseCells);
    Grid fineGrid(lowerLeft, upperRight, fineCells);

    auto coarseGridGeometry = std::make_shared<CoarseGridGeometry>(coarseGrid.leafGridView());
    auto fineGridGeometry = std::make_shared<GridGeometry>(fineGrid.leafGridView());
    const Dumux::TwoPVEColumnMapping<CoarseGridGeometry, GridGeometry, Scalar> mapping(coarseGridGeometry, fineGridGeometry);

    std::size_t expectedNumberOfColumns = 1;
    for (const auto cellsInDirection : coarseCells)
        expectedNumberOfColumns *= cellsInDirection;

    if (mapping.numberOfColumns() != expectedNumberOfColumns)
        DUNE_THROW(Dune::Exception, "Expected " << expectedNumberOfColumns << " columns, obtained " << mapping.numberOfColumns());

    for (std::size_t coarseIdx = 0; coarseIdx < mapping.numberOfColumns(); ++coarseIdx)
    {
        const auto& column = mapping.column(coarseIdx);
        if (column.size() != fineCells[dim-1])
            DUNE_THROW(Dune::Exception, "Expected " << fineCells[dim-1] << " fine elements in column " << coarseIdx << ", obtained " << column.size());

        for (std::size_t i = 0; i < column.size(); ++i)
        {
            const auto fineIdx = fineGridGeometry->elementMapper().index(column[i]);
            if (mapping.coarseIndex(fineIdx) != coarseIdx)
                DUNE_THROW(Dune::Exception, "Coarse-to-fine mapping is not invertible for fine element " << fineIdx);

            if (i > 0 && !(column[i-1].geometry().center()[dim-1] < column[i].geometry().center()[dim-1]))
                DUNE_THROW(Dune::Exception, "Column " << coarseIdx << " is not ordered vertically");
        }
    }

}

int main(int argc, char** argv)
{
    Dumux::initialize(argc, argv);
    checkMapping<2>();
    checkMapping<3>();
    return 0;
}
