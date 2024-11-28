// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test trace grid extraction from finite volume grid geometries.
 *
 */
#include <config.h>

#include <cmath>
#include <string_view>

#include <dune/common/exceptions.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>
#include <dumux/discretization/facetgrid.hh>


template<typename T>
std::vector<T> uniqueValuesIn(const std::vector<T>& in)
{
    std::vector<T> result = in;
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}


int main(int argc, char** argv)
{
    using namespace Dumux;
    initialize(argc, argv);

    using Grid = Dune::YaspGrid<2>;
    using GridView = typename Grid::LeafGridView;
    using GridGeometry = CCTpfaFVGridGeometry<GridView>;

    Grid grid{{1.0, 1.0}, {10, 10}};

    int exit_code = 0;
    const auto handleError = [&] (std::string_view message) {
        std::cout << message << std::endl;
        exit_code += 1;
    };

    {
        FVFacetGrid facetGrid{
            std::make_shared<GridGeometry>(grid.leafGridView()),
            [] (const auto& is) { return std::abs(is.geometry().center()[0] - 0.5) < 1e-6; }
        };
        if (facetGrid.gridView().size(0) != 10)
            handleError("Unexpected number of facet grid cells: " + std::to_string(facetGrid.gridView().size(0)));
        if (facetGrid.gridView().size(1) != 11)
            handleError("Unexpected number of facet grid vertices: " + std::to_string(facetGrid.gridView().size(1)));
    }

    {
        FVTraceGrid traceGrid{std::make_shared<GridGeometry>(grid.leafGridView())};
        if (traceGrid.gridView().size(0) != 40)
            handleError("Unexpected number of trace grid cells: " + std::to_string(traceGrid.gridView().size(0)));
        if (traceGrid.gridView().size(1) != 40)
            handleError("Unexpected number of trace grid vertices: " + std::to_string(traceGrid.gridView().size(0)));

        std::vector<std::size_t> adjacentElements;
        for (const auto& element: traceGrid.adjacentDomainElements())
            adjacentElements.push_back(grid.leafGridView().indexSet().index(element));
        if (adjacentElements.size() != 10 + 10 + 2*8)
            handleError("Unexpected number of adjacent domain elements");
        if (uniqueValuesIn(adjacentElements).size() != adjacentElements.size())
            handleError("Unexpected duplicates in adjacent domain elements");
    }

    return exit_code;
}
