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
#include <array>
#include <numeric>
#include <string_view>

#include <dune/common/float_cmp.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/alugrid/grid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>
#include <dumux/discretization/facetgrid.hh>


template<typename T>
std::vector<T> uniqueValuesIn(std::vector<T> in)
{
    std::sort(in.begin(), in.end());
    in.erase(std::unique(in.begin(), in.end()), in.end());
    return in;
}

template<typename FacetElementGeometry, typename DomainElementGeometry>
bool overlap(const FacetElementGeometry& facetGeo, const DomainElementGeometry& domainGeo)
{
    const auto cornersOf = [] (const auto& geo) {
        return std::views::iota(0, geo.corners())
            | std::views::transform([&] (int c) { return geo.corner(c); });
    };

    for (const auto& facetGeoCorner : cornersOf(facetGeo))
        if (std::ranges::none_of(cornersOf(domainGeo), [&] (const auto& domainGeoCorner) {
            return Dune::FloatCmp::eq(facetGeoCorner, domainGeoCorner);
        }))
            return false;

    return true;
}

template<typename FacetGridType>
int test()
{
    static constexpr int dim = FacetGridType::dimensionworld;
    static_assert(dim > 1);

    using Grid = Dune::YaspGrid<dim>;
    using GridView = typename Grid::LeafGridView;
    using GridGeometry = Dumux::CCTpfaFVGridGeometry<GridView>;

    const int cellsPerSide = 10;
    std::array<int, dim> cells; std::ranges::fill(cells, cellsPerSide);
    Dune::FieldVector<double, dim> size; std::ranges::fill(size, 1.0);

    Grid grid{size, cells};
    const auto gridView = grid.leafGridView();

    const auto cellsPerSlice = std::pow(cellsPerSide, dim-1);
    const auto pointPerSlice = std::pow(cellsPerSide+1, dim-1);
    const auto numBoundaryCells = cellsPerSlice*dim*2;
    const auto numBoundaryPoints = gridView.size(dim) - std::pow(cellsPerSide-1, dim);
    const auto numGridCellsTouchingBoundary = gridView.size(0) - std::pow(cellsPerSide-2, dim);

    int exitCode = 0;
    const auto handleError = [&] (std::string_view message) {
        std::cout << message <<  " @ dim = " << dim << std::endl;
        exitCode += 1;
    };

    {
        auto gridGeometry = std::make_shared<GridGeometry>(grid.leafGridView());
        auto facetGrid = makeFVFacetGrid<FacetGridType>(gridGeometry, [] (const auto& is) {
            return std::abs(is.geometry().center()[dim - 1] - 0.5) < 1e-6;
        });

        if (facetGrid.gridView().size(0) != cellsPerSlice)
            handleError("Unexpected number of facet grid cells: " + std::to_string(facetGrid.gridView().size(0)));
        if (facetGrid.gridView().size(dim-1) != pointPerSlice)
            handleError("Unexpected number of facet grid vertices: " + std::to_string(facetGrid.gridView().size(dim-1)));

        for (const auto& facetElement : elements(facetGrid.gridView()))
        {
            unsigned int count = 0;
            for (const auto& domainElement : facetGrid.domainElementsAdjacentTo(facetElement))
            {
                count++;
                if (!overlap(facetElement.geometry(), domainElement.geometry()))
                    handleError("Facet and domain element do not overlap");

                unsigned int scvfCount = 0;
                const auto fvGeometry = localView(*gridGeometry).bindElement(domainElement);
                for (const auto scvfIdx : facetGrid.domainScvfsAdjacentTo(facetElement, domainElement))
                {
                    scvfCount++;
                    const auto& scvf = fvGeometry.scvf(scvfIdx);
                    if (!overlap(fvGeometry.geometry(scvf), facetElement.geometry()))
                        handleError("Facet element and domain scvf do not overlap");
                    if (Dune::FloatCmp::ne(scvf.area(), facetElement.geometry().volume()))
                        handleError("Facet element and scvf areas do not match");
                }

                if (scvfCount != 1)
                    handleError("Expected one adjacent domain scvf per facet element");
            }

            if (count != 2)
                handleError("Expected two adjacent domain elements per facet element");
        }
    }

    {
        auto boundaryGrid = makeFVBoundaryGrid<FacetGridType>(std::make_shared<GridGeometry>(grid.leafGridView()));
        if (boundaryGrid.gridView().size(0) != numBoundaryCells)
            handleError("Unexpected number of trace grid cells: " + std::to_string(boundaryGrid.gridView().size(0)));
        if (boundaryGrid.gridView().size(dim-1) != numBoundaryPoints)
            handleError("Unexpected number of trace grid vertices: " + std::to_string(boundaryGrid.gridView().size(dim-1)));

        std::vector<std::size_t> adjacentElements;
        for (const auto& facetElement : elements(boundaryGrid.gridView()))
            for (const auto& element: boundaryGrid.domainElementsAdjacentTo(facetElement))
                adjacentElements.push_back(grid.leafGridView().indexSet().index(element));
        adjacentElements = uniqueValuesIn(adjacentElements);
        if (adjacentElements.size() != numGridCellsTouchingBoundary)
            handleError("Unexpected number of adjacent domain elements");
    }

    return exitCode;
}


int main(int argc, char** argv)
{
    using namespace Dumux;
    initialize(argc, argv);
    auto exitCode = test<Dune::FoamGrid<1, 2>>();
    exitCode += test<Dune::ALUGrid<2, 3, Dune::cube, Dune::nonconforming>>();
    return exitCode;
}
