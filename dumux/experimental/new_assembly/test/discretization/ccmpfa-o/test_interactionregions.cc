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
 * \brief Tests for the `CCMpfaO::InteractionRegions` class.
 */
#include <set>
#include <iostream>

#include <dune/common/timer.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/experimental/new_assembly/dumux/discretization/ccmpfa-o/interactionregions.hh>


template<typename Action>
auto printTime(const Action& action)
{
    Dune::Timer timer;
    auto result = action();
    timer.stop();
    std::cout << "Took " << timer.elapsed() << " seconds." << std::endl;
    return result;
}

template<typename Grid>
void test(const Grid& grid, const std::vector<std::size_t>& expectedSizes)
{
    std::cout << "\nConstructing interaction regions for grid with "
              << grid.leafGridView().size(0)
              << " cells in "
              << int(Grid::dimension) << "d" << std::endl;
    auto regions = printTime([&] () {
        return Dumux::CCMpfaO::InteractionRegions{
            grid.leafGridView(),
            grid.leafGridView().indexSet(),
            grid.leafGridView().indexSet()
        };
    });

    std::set<int> irSizes;
    for (std::size_t i = 0; i < regions.size(); ++i)
        irSizes.insert(regions[i].size());

    bool unexpected = std::ranges::any_of(irSizes, [&] (auto size) {
        return !std::ranges::count(expectedSizes, size);
    });
    if (unexpected)
        DUNE_THROW(Dune::InvalidStateException, "Found unexpected region size");
}

int main (int argc, char *argv[])
{
    test(Dune::YaspGrid<2>{{1.0, 1.0}, {20, 20}}, {1, 2, 4});
    test(Dune::YaspGrid<2>{{1.0, 1.0}, {200, 200}}, {1, 2, 4});
    test(Dune::YaspGrid<2>{{1.0, 1.0}, {2000, 2000}}, {1, 2, 4});

    test(Dune::YaspGrid<3>{{1.0, 1.0, 1.0}, {20, 20, 20}}, {1, 2, 4, 8});
    test(Dune::YaspGrid<3>{{1.0, 1.0, 1.0}, {200, 200, 200}}, {1, 2, 4, 8});

    std::cout << "\nAll tests passed" << std::endl;
    return 0;
}
