// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test facet grid extraction.
 *
 */
#include <cmath>
#include <array>
#include <string_view>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/alugrid/grid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/io/grid/facetgridmanager.hh>

template<typename T>
void removeDuplicates(std::vector<T>& in)
{
    std::sort(in.begin(), in.end());
    in.erase(std::unique(in.begin(), in.end()), in.end());
}

template<typename GridView>
void write(const GridView& gridView, const std::string& filename)
{
    Dune::VTKWriter{gridView}.write(filename);
}

template<typename FacetGrid>
int test()
{
    static constexpr int dim = FacetGrid::dimensionworld;
    static_assert(dim > 1);

    const int cellsPerSide = 10;
    std::array<int, dim> cells; std::ranges::fill(cells, cellsPerSide);
    Dune::FieldVector<double, dim> size; std::ranges::fill(size, 1.0);

    using Grid = Dune::YaspGrid<dim>;
    Grid grid{size, cells};
    write(grid.leafGridView(), "test_facetgrid_host_" + std::to_string(dim) + "d");

    const auto cellsPerSlice = std::pow(cellsPerSide, dim-1);
    const auto pointPerSlice = std::pow(cellsPerSide+1, dim-1);
    const auto numBoundaryCells = cellsPerSlice*dim*2;
    const auto numBoundaryPoints = grid.leafGridView().size(dim) - std::pow(cellsPerSide-1, dim);

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
        write(facetGridView, "test_facetgrid_interior_" + std::to_string(dim) + "d");

        if (facetGridView.size(0) != cellsPerSlice)
            handleError("Unexpected number of facet grid cells: " + std::to_string(facetGridView.size(0)));
        if (facetGridView.size(dim-1) != pointPerSlice)
            handleError("Unexpected number of facet grid vertices: " + std::to_string(facetGridView.size(dim-1)));

        std::vector<std::size_t> mappedHostGridVertexIndices;
        for (const auto& v : vertices(facetGridView))
            mappedHostGridVertexIndices.push_back(grid.leafGridView().indexSet().index(facetGridManager.hostGridVertex(v)));
        removeDuplicates(mappedHostGridVertexIndices);
        if (mappedHostGridVertexIndices.size() != pointPerSlice)
            handleError("Mapping from facet to host vertices is not unique");
    }

    { // grid composed of boundary facets
        Dumux::FacetGridManager<Grid, FacetGrid> facetGridManager;
        facetGridManager.init(grid, [] (const auto&, const auto& is) { return is.boundary(); });
        const auto& facetGridView = facetGridManager.grid().leafGridView();
        write(facetGridView, "test_facetgrid_boundary_" + std::to_string(dim) + "d");

        if (facetGridView.size(0) != numBoundaryCells)
            handleError("Unexpected number of trace grid cells: " + std::to_string(facetGridView.size(0)));
        if (facetGridView.size(dim-1) != numBoundaryPoints)
            handleError("Unexpected number of trace grid vertices: " + std::to_string(facetGridView.size(dim-1)));
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
