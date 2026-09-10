// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test facet grid extraction & the mapping of facet grid elements to the intersections,
 *        boundary faces, sub-control volume faces and degrees of freedom of a discretization.
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
#include <type_traits>

#include <dune/common/float_cmp.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/alugrid/grid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/geometry/geometryintersection.hh>
#include <dumux/geometry/intersectspointgeometry.hh>
#include <dumux/geometry/volume.hh>

#include <dumux/io/format.hh>
#include <dumux/io/grid/facetgridmanager.hh>
#include <dumux/discretization/facetgridmapper.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>
#include <dumux/discretization/pq1bubble/fvgridgeometry.hh>
#include <dumux/discretization/pq2/fvgridgeometry.hh>
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

template<typename Range>
auto toVector(Range&& range)
{
    std::vector<std::ranges::range_value_t<Range>> result;
    for (const auto& entry : range)
        result.push_back(entry);
    return result;
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
    using DM = typename GridGeometry::DiscretizationMethod;
    static constexpr bool isCVFE = Dumux::DiscretizationMethods::isCVFE<DM>;
    static constexpr bool isPQ2 = std::is_same_v<DM, Dumux::DiscretizationMethods::PQ2>;
    // the degrees of freedom on a facet are those of the Lagrange element of the same order on it
    const auto expectedNumDofs = isCVFE ? (isPQ2 ? std::pow(3, dim-1) : std::pow(2, dim-1)) : 1;

    const auto cellsPerSlice = std::pow(cellsPerSide, dim-1);
    const auto pointsPerSlice = std::pow(cellsPerSide+1, dim-1);
    const auto numBoundaryCells = cellsPerSlice*dim*2;
    const auto numBoundaryPoints = grid.leafGridView().size(dim) - std::pow(cellsPerSide-1, dim);
    const auto numGridCellsTouchingBoundary = grid.leafGridView().size(0) - std::pow(cellsPerSide-2, dim);

    int exitCode = 0;
    const auto handleError = [&] (std::string_view message) {
        std::cout << Dumux::Fmt::format("{} @ dim = {}", message, dim) << std::endl;
        exitCode += 1;
    };

    const auto testIntersections = [&] (const auto& facetElement,
                                        const auto& domainElement,
                                        const auto& mapper,
                                        const bool isBoundary = false) {
        const auto facetGeometry = facetElement.geometry();
        if (not intersectionVolume(facetGeometry, domainElement.geometry()).has_value())
            handleError("Facet and domain element do not overlap");

        const auto fvGeometry = localView(*gridGeometry).bindElement(domainElement);

        // the intersection the mapper reports is the one the facet element was extracted from
        {
            const auto isIdx = mapper.intersectionIndex(facetElement, domainElement);
            bool found = false;
            for (const auto& is : intersections(grid.leafGridView(), domainElement))
                if (is.indexInInside() == isIdx)
                {
                    found = true;
                    if ((is.geometry().center() - facetGeometry.center()).two_norm() > 1e-12)
                        handleError("The mapped intersection does not coincide with the facet element");
                }
            if (not found)
                handleError("The domain element has no intersection with the mapped index");
        }

        if (isBoundary)
        {
            const auto face = mapper.boundaryFace(fvGeometry, facetElement);
            if ((face.center() - facetGeometry.center()).two_norm() > 1e-12)
                handleError("The boundary face does not coincide with the facet element");
            if (Dune::FloatCmp::ne(face.area(), facetGeometry.volume(), 1e-7*facetGeometry.volume()))
                handleError("The boundary face does not have the area of the facet element");
        }

        const auto dofIndices = toVector(mapper.domainLocalDofsAdjacentTo(facetElement, domainElement));
        if (dofIndices.size() != expectedNumDofs)
            handleError(Dumux::Fmt::format(
                "Unexpected number of adjacent degrees of freedom: {}, expected {}",
                dofIndices.size(),
                expectedNumDofs
            ));

        if constexpr (isCVFE)
            for (const auto& localDof : localDofs(fvGeometry))
                if (std::ranges::find(dofIndices, localDof.index()) != dofIndices.end())
                    if (not Dumux::intersectsPointGeometry(ipData(fvGeometry, localDof).global(), facetGeometry))
                        handleError("An adjacent degree of freedom does not lie on the facet element");

        // the sub-control volume faces on a facet partition it; a control-volume finite element
        // scheme has none on an interior facet, since its faces lie inside the elements
        const auto scvfIndices = toVector(mapper.domainScvfsAdjacentTo(facetElement, domainElement));
        if (isCVFE and not isBoundary)
        {
            if (not scvfIndices.empty())
                handleError("Expected no sub-control volume faces on an interior facet");
        }
        else
        {
            double coveredArea = 0.0;
            for (const auto scvfIndex : scvfIndices)
            {
                const auto& scvf = fvGeometry.scvf(scvfIndex);
                const auto isVolume = intersectionVolume(facetGeometry, fvGeometry.geometry(scvf));
                if (not isVolume.has_value())
                    handleError("Facet element and domain scvf do not overlap");
                else
                    coveredArea += isVolume.value();
            }
            if (Dune::FloatCmp::ne(coveredArea, facetGeometry.volume(), 1e-7*facetGeometry.volume()))
                handleError(Dumux::Fmt::format(
                    "The sub-control volume faces cover {} of the facet element, expected {}",
                    coveredArea,
                    facetGeometry.volume()
                ));
        }
    };

    { // grid composed of interior facets
        Dumux::FacetGridManager<Grid, FacetGrid> facetGridManager;
        facetGridManager.init(grid, [] (const auto&, const auto& is) {
            return std::abs(is.geometry().center()[dim - 1] - 0.5) < 1e-6;
        });
        const auto& facetGridView = facetGridManager.grid().leafGridView();
        Dumux::FacetGridMapper mapper{facetGridManager, gridGeometry};

        if (facetGridView.size(0) != cellsPerSlice)
            handleError(Dumux::Fmt::format("Unexpected number of facet grid cells: {}", facetGridView.size(0)));
        if (facetGridView.size(dim-1) != pointsPerSlice)
            handleError(Dumux::Fmt::format("Unexpected number of facet grid vertices: {}", facetGridView.size(dim-1)));

        for (const auto& facetElement : elements(facetGridView))
        {
            unsigned int elementCount = 0;
            for (const auto& domainElement : mapper.domainElementsAdjacentTo(facetElement))
            {
                elementCount++;
                testIntersections(facetElement, domainElement, mapper);
            }

            if (elementCount != 2)
                handleError(Dumux::Fmt::format("Expected two adjacent domain elements per facet element, found {}", elementCount));
        }
    }

    { // grid composed of boundary facets
        Dumux::FacetGridManager<Grid, FacetGrid> facetGridManager;
        facetGridManager.init(grid, [] (const auto&, const auto& is) { return is.boundary(); });
        const auto& facetGridView = facetGridManager.grid().leafGridView();
        Dumux::FacetGridMapper mapper{facetGridManager, gridGeometry};

        if (facetGridView.size(0) != numBoundaryCells)
            handleError(Dumux::Fmt::format("Unexpected number of trace grid cells: {}", facetGridView.size(0)));
        if (facetGridView.size(dim-1) != numBoundaryPoints)
            handleError(Dumux::Fmt::format("Unexpected number of trace grid vertices: {}", facetGridView.size(dim-1)));

        std::vector<std::size_t> adjacentElements;
        for (const auto& facetElement : elements(facetGridView))
            for (const auto& element : mapper.domainElementsAdjacentTo(facetElement))
            {
                adjacentElements.push_back(grid.leafGridView().indexSet().index(element));
                testIntersections(facetElement, element, mapper, true);
            }

        std::ranges::sort(adjacentElements);
        adjacentElements.erase(std::unique(adjacentElements.begin(), adjacentElements.end()), adjacentElements.end());
        if (adjacentElements.size() != numGridCellsTouchingBoundary)
            handleError("Unexpected number of adjacent domain elements");
    }

    return exitCode;
}

template<typename GV> using TpfaGridGeometry = Dumux::CCTpfaFVGridGeometry<GV>;
template<typename GV> using BoxGridGeometry = Dumux::BoxFVGridGeometry<double, GV>;
template<typename GV> using PQ1BubbleGridGeometry = Dumux::PQ1BubbleFVGridGeometry<double, GV>;
template<typename GV> using PQ2GridGeometry = Dumux::PQ2FVGridGeometry<double, GV>;

int main(int argc, char** argv)
{
    using namespace Dumux;

    int exitCode = 0;
    initialize(argc, argv);
    {
        using Grid = Dune::FoamGrid<1, 2>;
        exitCode += test<Grid, TpfaGridGeometry>();
        exitCode += test<Grid, BoxGridGeometry>();
        exitCode += test<Grid, PQ1BubbleGridGeometry>();
        exitCode += test<Grid, PQ2GridGeometry>();
    }
    {
        using Grid = Dune::ALUGrid<2, 3, Dune::cube, Dune::nonconforming>;
        exitCode += test<Grid, TpfaGridGeometry>();
        exitCode += test<Grid, BoxGridGeometry>();
        exitCode += test<Grid, PQ1BubbleGridGeometry>();
        exitCode += test<Grid, PQ2GridGeometry>();
    }
    return exitCode;
}
