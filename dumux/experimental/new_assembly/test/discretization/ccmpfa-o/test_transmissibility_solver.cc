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
 * \brief Tests for the `CCMpfaO::TransmissibilitySolver` class.
 */
#include <iostream>
#include <vector>

#include <dune/common/timer.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/experimental/new_assembly/dumux/discretization/ccmpfa-o/interactionregions.hh>
#include <dumux/experimental/new_assembly/dumux/discretization/ccmpfa-o/transmissibilitysolver.hh>

template<typename GV>
class DummyGridGeometry
{
    using Element = typename GV::template Codim<0>::Entity;

public:
    using GridView = GV;

    explicit DummyGridGeometry(const GridView& gv)
    : gridView_(gv)
    {
        elementMap_.reserve(gv.size(0));
        for (auto element : elements(gv))
            elementMap_.emplace_back(std::move(element));
    }

    const GridView& gridView() const { return gridView_; }
    const Element& element(std::size_t i) const { return elementMap_[i]; }
    const auto& elementMapper() const { return gridView_.indexSet(); }
    const auto& vertexMapper() const { return gridView_.indexSet(); }

private:
    const GridView gridView_;
    std::vector<Element> elementMap_;
};


template<typename Action>
void printTime(const Action& action)
{
    Dune::Timer timer;
    action();
    timer.stop();
    std::cout << "Took " << timer.elapsed() << " seconds." << std::endl;
}

template<typename GridGeometry>
auto makeSolvers(const GridGeometry& gridGeometry)
{
    using Dumux::CCMpfaO::InteractionRegions;
    using Dumux::CCMpfaO::TransmissibilitySolver;

    InteractionRegions regions{gridGeometry.gridView(),
                               gridGeometry.elementMapper(),
                               gridGeometry.vertexMapper()};

    std::vector<TransmissibilitySolver<GridGeometry>> result;
    result.reserve(gridGeometry.gridView().size(GridGeometry::GridView::dimension));
    for (const auto& v : vertices(gridGeometry.gridView()))
        result.emplace_back(TransmissibilitySolver{
            gridGeometry,
            regions[gridGeometry.vertexMapper().index(v)]
        });

    return result;
}

template<typename Grid>
void test(const Grid& grid)
{
    using Dumux::CCMpfaO::InteractionRegions;
    using Dumux::CCMpfaO::TransmissibilitySolver;

    DummyGridGeometry gridGeometry{grid.leafGridView()};
    printTime([&] () { makeSolvers(gridGeometry); });
}

int main (int argc, char *argv[])
{
    const unsigned int numCells = 10;
    test(Dune::YaspGrid<2>{
        {1.0, 1.0},
        {numCells, numCells}
    });

    std::cout << "\nAll tests passed" << std::endl;
    return 0;
}
