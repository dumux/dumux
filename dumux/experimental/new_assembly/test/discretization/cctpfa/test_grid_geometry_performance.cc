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
// TODO: DELETE THIS TEST
#include <iostream>
#include <concepts>
#include <cmath>

#include <dune/common/timer.hh>
#include <dune/common/exceptions.hh>
#include <dune/grid/yaspgrid.hh>

#ifndef USE_OLD
#define USE_OLD 0
#endif

#ifndef USE_CACHING
#define USE_CACHING 0
#endif

#if USE_OLD
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>
#else
#include <dumux/experimental/new_assembly/dumux/discretization/cctpfa/gridgeometry.hh>
#endif

static constexpr int numFaceAccesses = 5;
static constexpr bool caching = USE_CACHING;

template<std::floating_point T1, std::floating_point T2>
bool equal(T1 t1, T2 t2)
{
    using std::abs;
    return abs(t1 - t2) < 1e-7;
}

template<typename Action>
decltype(auto) printTime(const std::string& name, const Action& action)
{
    static constexpr bool isVoidAction = std::same_as<
        void, std::invoke_result_t<Action>
    >;

    Dune::Timer timer;
    if constexpr (isVoidAction)
    {
        action();
        timer.stop();
        std::cout << name << " took " << timer.elapsed() << " seconds." << std::endl;
    }
    else
    {
        auto result = action();
        timer.stop();
        std::cout << name << " took " << timer.elapsed() << " seconds." << std::endl;
        return result;
    }
}

template<template<typename GV> typename GridGeometry,
         typename Helper,
         typename Grid>
void runTests(const Grid& grid, std::size_t numCellsPerDirection)
{
    const auto& gridView = grid.leafGridView();
    const double expectedScvfArea = 1.0/numCellsPerDirection;
    const double expectedScvVolume = expectedScvfArea*expectedScvfArea;
    const std::size_t numExpectedScvfsPerScv = Grid::dimension == 2 ? 4 : 6;
    std::size_t expectedNumScvf = 0;
    for (const auto& e : elements(gridView))
        expectedNumScvf += e.subEntities(1);

    std::cout << "Constructing grid geometry... ";
    auto gridGeometry = printTime("GridGeometry construction", [&] () {
        return GridGeometry<typename Grid::LeafGridView>{gridView};
    });

    if (gridGeometry.numDofs() != gridView.size(0))
        DUNE_THROW(Dune::InvalidStateException, "Unexpected number of dofs");

    printTime("Running the tests", [&] () {
        auto fvGeometry = localView(gridGeometry);
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            fvGeometry.bind(element);
            std::size_t insideDofIndex = std::numeric_limits<std::size_t>::max();
            for (const auto& scv : scvs(fvGeometry))
            {
                if (!equal(scv.volume(), expectedScvVolume))
                    DUNE_THROW(Dune::InvalidStateException, "Unexpected scv volume");
                insideDofIndex = scv.dofIndex();
            }

            int scvfCount = 0;
            for (const auto& scvf : scvfs(fvGeometry))
            {
                scvfCount++;
                if (!equal(scvf.area(), expectedScvfArea))
                    DUNE_THROW(Dune::InvalidStateException, "Unexpected scvf area");
                if (Helper::insideScv(fvGeometry, scvf).dofIndex() != insideDofIndex)
                    DUNE_THROW(Dune::InvalidStateException, "Unexpected inside scv");
                if (!Helper::onBoundary(fvGeometry, scvf))
                {
                    // get the normal vector multiple times to mimmick local assembly
                    for (int i = 0; i < numFaceAccesses; ++i)
                    {
                        const auto& insideScv = Helper::insideScv(fvGeometry, scvf);
                        const auto& outsideScv = Helper::outsideScv(fvGeometry, scvf);
                        const auto n = scvf.unitOuterNormal();
                        auto d = outsideScv.center() - insideScv.center();
                        d /= d.two_norm();
                        if ((d-n).two_norm() > 1e-7)
                            DUNE_THROW(Dune::InvalidStateException, "Unexpected normal/outsideScv alignment");
                    }
                }
                else
                {
                    // get the normal vector multiple times to mimmick local assembly
                    for (int i = 0; i < numFaceAccesses; ++i)
                    {
                        const auto n = scvf.unitOuterNormal();
                        auto d = scvf.center() -  Helper::insideScv(fvGeometry, scvf).center();
                        d /= d.two_norm();
                        if ((d-n).two_norm() > 1e-7)
                            DUNE_THROW(Dune::InvalidStateException, "Unexpected scvf/cell center alignment");
                    }
                }
            }
            if (scvfCount != numExpectedScvfsPerScv)
                DUNE_THROW(Dune::InvalidStateException, "Unexpected number of scvfs per scv");
        }
    });
}

#if USE_OLD
template<typename GV>
using GridGeometry = Dumux::CCTpfaFVGridGeometry<GV, caching>;

struct GridGeometryHelper
{
    template<typename LocalView, typename Scvf>
    static const auto& insideScv(const LocalView& lv, const Scvf& scvf)
    { return lv.scv(scvf.insideScvIdx()); }

    template<typename LocalView, typename Scvf>
    static const auto& outsideScv(const LocalView& lv, const Scvf& scvf)
    { return lv.scv(scvf.outsideScvIdx()); }

    template<typename LocalView, typename Scvf>
    static bool onBoundary(const LocalView& lv, const Scvf& scvf)
    { return scvf.boundary(); }
};

std::string prefix() { return "old"; }
#else
template<typename GV>
using GridGeometry = Dumux::CCTpfaGridGeometry<GV, caching>;

struct GridGeometryHelper
{
    template<typename LocalView, typename Scvf>
    static const auto& insideScv(const LocalView& lv, const Scvf& scvf)
    { return lv.insideScv(scvf); }

    template<typename LocalView, typename Scvf>
    static const auto& outsideScv(const LocalView& lv, const Scvf& scvf)
    { return lv.outsideScv(scvf); }

    template<typename LocalView, typename Scvf>
    static bool onBoundary(const LocalView& lv, const Scvf& scvf)
    { return lv.onBoundary(scvf); }
};

std::string prefix() { return "new"; }
#endif

#if USE_CACHING
std::string suffix() { return "with caching"; }
#else
std::string suffix() { return "without caching"; }
#endif

int main (int argc, char *argv[])
{
    const unsigned int numCells = 2500;
    Dune::YaspGrid<2> grid{
        {1.0, 1.0},
        {numCells, numCells}
    };

    printTime("Testing " + prefix() + " grid geometry " + suffix(), [&] () {
        runTests<GridGeometry, GridGeometryHelper>(grid, numCells);
    });

    std::cout << "\nAll tests passed" << std::endl;
    return 0;
}
