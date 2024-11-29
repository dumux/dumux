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

#include <dune/common/float_cmp.hh>
#include <dune/grid/yaspgrid.hh>

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


int main(int argc, char** argv)
{
    using namespace Dumux;
    initialize(argc, argv);

    using Grid = Dune::YaspGrid<2>;
    using GridView = typename Grid::LeafGridView;
    using GridGeometry = CCTpfaFVGridGeometry<GridView>;

    Grid grid{{1.0, 1.0}, {10, 10}};

    int exitCode = 0;
    const auto handleError = [&] (std::string_view message) {
        std::cout << message << std::endl;
        exitCode += 1;
    };

    {
        auto gridGeometry = std::make_shared<GridGeometry>(grid.leafGridView());
        FVFacetGrid facetGrid{gridGeometry, [] (const auto& is) {
            return std::abs(is.geometry().center()[0] - 0.5) < 1e-6;
        }};

        if (facetGrid.gridView().size(0) != 10)
            handleError("Unexpected number of facet grid cells: " + std::to_string(facetGrid.gridView().size(0)));
        if (facetGrid.gridView().size(1) != 11)
            handleError("Unexpected number of facet grid vertices: " + std::to_string(facetGrid.gridView().size(1)));

        for (const auto& facetElement : elements(facetGrid.gridView()))
        {
            unsigned int count = 0;
            for (const auto& domainElement : facetGrid.domainElementsAdjacentTo(facetElement))
            {
                count++;
                if (!overlap(facetElement.geometry(), domainElement.geometry()))
                    handleError("Facet and Domain element do not overlap");

                unsigned int scvfCount = 0;
                const auto fvGeometry = localView(*gridGeometry).bindElement(domainElement);
                for (const auto scvfIdx : facetGrid.domainScvfsAdjacentTo(facetElement, domainElement))
                {
                    scvfCount++;
                    if (!overlap(fvGeometry.geometry(fvGeometry.scvf(scvfIdx)), facetElement.geometry()))
                        handleError("Facet element and domain scvf do not overlap");
                }

                if (scvfCount != 1)
                    handleError("Expected one adjacent domain scvf per facet element");
            }

            if (count != 2)
                handleError("Expected two adjacent domain elements per facet element");
        }
    }

    {
        FVTraceGrid traceGrid{std::make_shared<GridGeometry>(grid.leafGridView())};
        if (traceGrid.gridView().size(0) != 40)
            handleError("Unexpected number of trace grid cells: " + std::to_string(traceGrid.gridView().size(0)));
        if (traceGrid.gridView().size(1) != 40)
            handleError("Unexpected number of trace grid vertices: " + std::to_string(traceGrid.gridView().size(0)));

        std::vector<std::size_t> adjacentElements;
        for (const auto& facetElement : elements(traceGrid.gridView()))
            for (const auto& element: traceGrid.domainElementsAdjacentTo(facetElement))
                adjacentElements.push_back(grid.leafGridView().indexSet().index(element));
        adjacentElements = uniqueValuesIn(adjacentElements);
        if (adjacentElements.size() != 10 + 10 + 2*8)
            handleError("Unexpected number of adjacent domain elements");
    }

    return exitCode;
}
