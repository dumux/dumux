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
#include <vector>
#include <ranges>
#include <numeric>
#include <algorithm>
#include <string_view>
#include <iterator>
#include <optional>

#include <dune/common/float_cmp.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/alugrid/grid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/geometry/geometryintersection.hh>
#include <dumux/geometry/volume.hh>

#include <dumux/io/grid/facetgridmanager.hh>
#include <dumux/discretization/facetgridmapper.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>


template<typename FacetElementGeometry, typename DomainElementGeometry>
std::optional<double> intersectionVolume(const FacetElementGeometry& facetGeo, const DomainElementGeometry& domainGeo)
{
    using Algorithm = Dumux::GeometryIntersection<FacetElementGeometry, DomainElementGeometry>;
    typename Algorithm::Intersection intersection;
    if (not Algorithm::intersection(facetGeo, domainGeo, intersection))
        return {};

    static constexpr int facetDim = FacetElementGeometry::mydimension;
    static_assert(facetDim == 1 or facetDim == 2);
    if constexpr (facetDim == 1)
        return (intersection[1] - intersection[0]).two_norm();
    else
    {
        if (intersection.size() != 4)
            DUNE_THROW(Dune::InvalidStateException, "Expected quadrilateral intersection");
        return Dumux::convexPolytopeVolume<facetDim>(
            Dune::GeometryTypes::quadrilateral,
            [&] (unsigned int i) {
                static constexpr int map[4] = {0, 1, 3, 2};
                return intersection[map[i]];
            }
        );
    }
}

template<typename FacetGrid, template<typename> typename GG>
int test()
{
    static constexpr int dim = FacetGrid::dimensionworld;
    static_assert(dim > 1);
    static_assert(dim == int(FacetGrid::dimension) + 1);

    const int cellsPerSide = 10;
    std::array<int, dim> cells; std::ranges::fill(cells, cellsPerSide);
    Dune::FieldVector<double, dim> size; std::ranges::fill(size, 1.0);

    using Grid = Dune::YaspGrid<dim>;
    using GridGeometry = GG<typename Grid::LeafGridView>;
    Grid grid{size, cells};
    auto gridGeometry = std::make_shared<GridGeometry>(grid.leafGridView());
    static constexpr bool isBox = Dumux::DiscretizationMethods::isCVFE<typename GridGeometry::DiscretizationMethod>;

    const auto cellsPerSlice = std::pow(cellsPerSide, dim-1);
    const auto pointsPerSlice = std::pow(cellsPerSide+1, dim-1);
    const auto numBoundaryCells = cellsPerSlice*dim*2;
    const auto numBoundaryPoints = grid.leafGridView().size(dim) - std::pow(cellsPerSide-1, dim);
    const auto numGridCellsTouchingBoundary = grid.leafGridView().size(0) - std::pow(cellsPerSide-2, dim);

    int exitCode = 0;
    const auto handleError = [&] (std::string_view message) {
        std::cout << message <<  " @ dim = " << dim << std::endl;
        exitCode += 1;
    };

    { // grid composed of interior facets
        Dumux::FacetGridManager<Grid, FacetGrid> facetGridManager;
        facetGridManager.init(grid, [] (const auto&, const auto& is) {
            return std::abs(is.geometry().center()[dim - 1] - 0.5) < 1e-6;
        });
        const auto& facetGridView = facetGridManager.grid().leafGridView();
        Dumux::FVFacetGridMapper mapper{facetGridView, gridGeometry};

        if (facetGridView.size(0) != cellsPerSlice)
            handleError("Unexpected number of facet grid cells: " + std::to_string(facetGridView.size(0)));
        if (facetGridView.size(dim-1) != pointsPerSlice)
            handleError("Unexpected number of facet grid vertices: " + std::to_string(facetGridView.size(dim-1)));

        for (const auto& facetElement : elements(facetGridView))
        {
            unsigned int elementCount = 0;
            for (const auto& domainElement : mapper.domainElementsAdjacentTo(facetElement))
            {
                elementCount++;
                if (not intersectionVolume(facetElement.geometry(), domainElement.geometry()).has_value())
                    handleError("Facet and domain element do not overlap");

                const auto sceIndices = [&] () {
                    std::vector<std::size_t> result;
                    if constexpr (isBox)
                        std::ranges::copy(mapper.domainScvsAdjacentTo(facetElement, domainElement), std::back_inserter(result));
                    else
                        std::ranges::copy(mapper.domainScvfsAdjacentTo(facetElement, domainElement), std::back_inserter(result));
                    return result;
                } ();

                const auto expectedNumSces = isBox ? std::pow(2, dim-1) : 1;
                if (sceIndices.size() != expectedNumSces)
                    handleError(
                        "Unexpected number of adjacent sub-control entities: " + std::to_string(sceIndices.size())
                        + " expected " + std::to_string(expectedNumSces)
                    );

                const auto fvGeometry = localView(*gridGeometry).bindElement(domainElement);
                for (const auto sceIndex : sceIndices)
                {
                    const auto sce = [&] () {
                        if constexpr (isBox) return fvGeometry.scv(sceIndex);
                        else return fvGeometry.scvf(sceIndex);
                    } ();
                    const auto isVolume = intersectionVolume(facetElement.geometry(), fvGeometry.geometry(sce));
                    const auto expectedVolume = facetElement.geometry().volume()/expectedNumSces;
                    if (not isVolume.has_value())
                        handleError("Facet element and domain sce do not overlap");
                    if (isVolume.has_value() and Dune::FloatCmp::ne(isVolume.value(), expectedVolume, expectedVolume*1e-7))
                        handleError(
                            "Unexpected area/volume of the intersection between sub-control entity and facet element: "
                            + std::to_string(isVolume.value()) + " vs. " + std::to_string(expectedVolume)
                        );
                }
            }

            if (elementCount != 2)
                handleError("Expected two adjacent domain elements per facet element, found " + std::to_string(elementCount));
        }
    }

    { // grid composed of boundary facets
        Dumux::FacetGridManager<Grid, FacetGrid> facetGridManager;
        facetGridManager.init(grid, [] (const auto&, const auto& is) { return is.boundary(); });
        const auto& facetGridView = facetGridManager.grid().leafGridView();
        Dumux::FVFacetGridMapper mapper{facetGridView, gridGeometry};

        if (facetGridView.size(0) != numBoundaryCells)
            handleError("Unexpected number of trace grid cells: " + std::to_string(facetGridView.size(0)));
        if (facetGridView.size(dim-1) != numBoundaryPoints)
            handleError("Unexpected number of trace grid vertices: " + std::to_string(facetGridView.size(dim-1)));

        std::vector<std::size_t> adjacentElements;
        for (const auto& facetElement : elements(facetGridView))
            for (const auto& element : mapper.domainElementsAdjacentTo(facetElement))
                adjacentElements.push_back(grid.leafGridView().indexSet().index(element));

        std::ranges::sort(adjacentElements);
        adjacentElements.erase(std::unique(adjacentElements.begin(), adjacentElements.end()), adjacentElements.end());
        if (adjacentElements.size() != numGridCellsTouchingBoundary)
            handleError("Unexpected number of adjacent domain elements");
    }

    return exitCode;
}

template<typename GV> using TpfaGridGeometry = Dumux::CCTpfaFVGridGeometry<GV>;
template<typename GV> using BoxGridGeometry = Dumux::BoxFVGridGeometry<double, GV>;

int main(int argc, char** argv)
{
    using namespace Dumux;

    int exitCode = 0;
    initialize(argc, argv);
    {
        using Grid = Dune::FoamGrid<1, 2>;
        exitCode += test<Grid, TpfaGridGeometry>();
        exitCode += test<Grid, BoxGridGeometry>();
    }
    {
        using Grid = Dune::ALUGrid<2, 3, Dune::cube, Dune::nonconforming>;
        exitCode += test<Grid, TpfaGridGeometry>();
        exitCode += test<Grid, BoxGridGeometry>();
    }
    return exitCode;
}
