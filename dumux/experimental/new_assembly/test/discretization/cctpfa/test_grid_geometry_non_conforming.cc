// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup CCTpfaDiscretization
 * \brief Test for the `CCTpfaGridGeometry` class on a locally refined grid.
 */
#include <config.h>

#include <iostream>
#include <algorithm>
#include <cmath>

#include <dune/alugrid/grid.hh>
#include <dumux/common/initialize.hh>

#include <dumux/experimental/new_assembly/dumux/discretization/cctpfa/gridgeometry.hh>

static constexpr double domainSize = 1.0;
static constexpr int numCellsPerSide = 3;

template<typename Element>
bool isCenterElement(const Element& element)
{
    using std::abs;
    const auto center = element.geometry().center();
    return abs(center[0] - domainSize/2.0) < 1e-6
        && abs(center[1] - domainSize/2.0) < 1e-6;
}

template<typename Grid>
void refineCenterElement(Grid& grid)
{
    for (const auto& element : elements(grid.leafGridView()))
        if (isCenterElement(element))
            grid.mark(1, element);

    grid.preAdapt();
    if (!grid.adapt())
        DUNE_THROW(Dune::InvalidStateException, "Adaptation did not occur");
    grid.postAdapt();
}

template<typename Grid>
auto makeGrid()
{
    std::unique_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(
        {0.0, 0.0}, {domainSize,domainSize}, {numCellsPerSide, numCellsPerSide}
    );

    refineCenterElement(*grid);
    if (grid->leafGridView().size(0) != 12)
        DUNE_THROW(Dune::InvalidStateException, "Unexpected number of cells after adaptation: " << grid->leafGridView().size(0));

    return grid;
}

template<typename GridView, typename Element>
std::vector<double> expectedScvfAreas(const GridView& gridView, const Element& element)
{
    std::vector<double> areas;
    for (const auto& is : intersections(gridView, element))
        areas.push_back(is.geometry().volume());
    return areas;
}

template<typename GridView, typename Element>
auto expectedScvfCenters(const GridView& gridView, const Element& element)
{
    std::vector<typename Element::Geometry::GlobalCoordinate> centers;
    for (const auto& is : intersections(gridView, element))
        centers.push_back(is.geometry().center());
    return centers;
}

template<typename GridGeometry, typename Element>
void checkScvfAreas(const GridGeometry& gridGeometry, const Element& element)
{
    const auto fvGeometry = localView(gridGeometry).bindElement(element);
    auto expectedAreas = expectedScvfAreas(gridGeometry.gridView(), element);
    std::ranges::sort(expectedAreas);

    const auto testScvfAreas = [&] (auto& scvfAreas) {
        if (expectedAreas.size() != scvfAreas.size())
            DUNE_THROW(Dune::InvalidStateException, "Unexpected number of scvfs");

        using std::abs;
        std::ranges::sort(scvfAreas);
        for (std::size_t i = 0; i < scvfAreas.size(); ++i)
            if (abs(scvfAreas[i] - expectedAreas[i]) > 1e-6)
                DUNE_THROW(Dune::InvalidStateException, "Unexpected scvf area");
    };

    std::vector<double> scvfAreas;
    for (const auto& scvf : scvfs(fvGeometry))
        scvfAreas.push_back(scvf.area());
    testScvfAreas(scvfAreas);

    scvfAreas.clear();
    for (const auto& scvf : scvfs(fvGeometry))
        scvfAreas.push_back(fvGeometry.geometry(scvf).volume());
    testScvfAreas(scvfAreas);
}

template<typename GridGeometry, typename Element>
void checkScvfCenters(const GridGeometry& gridGeometry, const Element& element)
{
    auto expectedCenters = expectedScvfCenters(gridGeometry.gridView(), element);

    std::vector<typename Element::Geometry::GlobalCoordinate> scvfCenters;
    for (const auto& scvf : scvfs(localView(gridGeometry).bind(element)))
        scvfCenters.push_back(scvf.center());

    std::ranges::all_of(scvfCenters, [&] (const auto& center) {
        return std::ranges::any_of(expectedCenters, [&] (const auto& expectedCenter) {
            return (expectedCenter - center).two_norm() < 1e-6;
        });
    });
}


template<typename GridGeometry, typename Element>
void checkElementScvfs(const GridGeometry& gridGeometry, const Element& element)
{
    checkScvfAreas(gridGeometry, element);
    checkScvfCenters(gridGeometry, element);
}


int main (int argc, char *argv[])
{
    Dumux::initialize(argc, argv);

    auto grid = makeGrid<Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>>();
    Dumux::CCTpfaGridGeometry gridGeometry{grid->leafGridView()};
    if (gridGeometry.numScv() != 12)
        DUNE_THROW(Dune::InvalidStateException, "Unexpected number of scvs: " << gridGeometry.numScv());
    if (gridGeometry.numScvf() != 52)
        DUNE_THROW(Dune::InvalidStateException, "Unexpected number of scvfs: " << gridGeometry.numScvf());

    std::vector<std::size_t> scvfCounts;
    for (const auto& element : elements(gridGeometry.gridView()))
        checkElementScvfs(gridGeometry, element);

    std::cout << "\nAll tests passed" << std::endl;
    return 0;
}
