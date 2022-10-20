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
 * \brief Test for the `CCTpfaGridGeometry` class on a network grid.
 */
#include <config.h>

#include <iostream>
#include <algorithm>
#include <cmath>

#include <dune/foamgrid/foamgrid.hh>
#include <dumux/common/initialize.hh>

#include <dumux/experimental/new_assembly/dumux/discretization/cctpfa/gridgeometry.hh>

template<typename Grid>
auto makeGrid()
{
    Dune::GridFactory<Grid> factory;
    factory.insertVertex({-1.0, -1.0});
    factory.insertVertex({0.0, 0.0});
    factory.insertVertex({1.0, -1.0});
    factory.insertVertex({0.0, 1.0});
    factory.insertElement(Dune::GeometryTypes::line, {0, 1});
    factory.insertElement(Dune::GeometryTypes::line, {1, 2});
    factory.insertElement(Dune::GeometryTypes::line, {1, 3});
    return factory.createGrid();
}

template<typename GlobalPosition>
bool isBifurcation(const GlobalPosition& pos)
{
    return std::all_of(pos.begin(), pos.end(), [] (auto x) {
        using std::abs;
        return abs(x) < 1e-6;
    });
}

template<typename GridGeometry>
void testGridGeometry(const GridGeometry& gridGeometry)
{
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        auto fvGeometry = localView(gridGeometry).bind(element);
        for (const auto& scvf : scvfs(fvGeometry))
            if (isBifurcation(scvf.ipGlobal()))
                if (fvGeometry.numOutsideNeighbors(scvf) != 2) {
                    std::cout<<"asid"<<std::endl;
                    for (int i = 0; i < fvGeometry.numOutsideNeighbors(scvf); ++i)
                        std::cout << "O = " << fvGeometry.outsideScv(scvf, i).dofIndex() << std::endl;
                    DUNE_THROW(Dune::InvalidStateException, "Unexpected number of outside Scvs: " << fvGeometry.numOutsideNeighbors(scvf));
                }
    }
}

int main (int argc, char *argv[])
{
    Dumux::initialize(argc, argv);

    using Grid = Dune::FoamGrid<1, 2>;
    using GridView = typename Grid::LeafGridView;

    auto grid = makeGrid<Grid>();
    if (grid->leafGridView().size(0) != 3)
        DUNE_THROW(Dune::InvalidStateException, "Unexpected number of elements:" << grid->leafGridView().size(0));
    if (grid->leafGridView().size(1) != 4)
        DUNE_THROW(Dune::InvalidStateException, "Unexpected number of vertices: " << grid->leafGridView().size(1));

    testGridGeometry(Dumux::CCTpfaGridGeometry<GridView, true>{grid->leafGridView()});
    testGridGeometry(Dumux::CCTpfaGridGeometry<GridView, false>{grid->leafGridView()});

    std::cout << "\nAll tests passed" << std::endl;
    return 0;
}
