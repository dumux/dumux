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
 * \brief Tests for the `CCMpfaOGridGeometry` class.
 */
#include <iostream>
#include <concepts>
#include <cmath>
#include <thread>

#include <dune/common/timer.hh>
#include <dune/common/exceptions.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>
#include <dumux/experimental/new_assembly/dumux/discretization/ccmpfa-o/gridgeometry.hh>

static constexpr int numFaceAccesses = 1;

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

template<typename Grid>
void testStructuredGrid2D(const Grid& grid, std::size_t numCellsPerDirection)
{
    const auto& gridView = grid.leafGridView();
    const double expectedScvfArea = 1.0/numCellsPerDirection;
    const double expectedScvVolume = expectedScvfArea*expectedScvfArea;
    const std::size_t numExpectedScvfsPerScv = 4;
    std::size_t expectedNumScvf = 0;
    for (const auto& e : elements(gridView))
        expectedNumScvf += e.subEntities(1);

    std::cout << "Constructing grid geometry... ";
    auto gridGeometry = printTime("GridGeometry construction", [&] () {
        return Dumux::CCMpfaOGridGeometry{gridView};
    });

    if (gridGeometry.numDofs() != gridView.size(0))
        DUNE_THROW(Dune::InvalidStateException, "Unexpected number of dofs");
    if (gridGeometry.numScv() != gridView.size(0))
        DUNE_THROW(Dune::InvalidStateException,
                   "Unexpected number of scvs " << gridGeometry.numScv() << " (expected " << gridView.size(0) << ")");
    if (gridGeometry.numScvf() != expectedNumScvf)
        DUNE_THROW(Dune::InvalidStateException,
                   "Unexpected number of scvfs " << gridGeometry.numScv() << " (expected " << expectedNumScvf << ")");

    printTime("Running the tests", [&] () {
        auto localGridGeometry = localView(gridGeometry);
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            localGridGeometry.bind(element);
            std::size_t insideDofIndex = std::numeric_limits<std::size_t>::max();
            for (const auto& scv : scvs(localGridGeometry))
            {
                if (!equal(localGridGeometry.volume(scv), expectedScvVolume))
                    DUNE_THROW(Dune::InvalidStateException, "Unexpected scv volume");
                if (!equal(localGridGeometry.geometry(scv).volume(), expectedScvVolume))
                    DUNE_THROW(Dune::InvalidStateException, "Unexpected scv geometry volume");
                insideDofIndex = localGridGeometry.dofIndex(scv);
            }

            int scvfCount = 0;
            for (const auto& scvf : scvfs(localGridGeometry))
            {
                scvfCount++;
                if (!equal(localGridGeometry.area(scvf), expectedScvfArea))
                    DUNE_THROW(Dune::InvalidStateException, "Unexpected scvf area");
                if (!equal(localGridGeometry.geometry(scvf).volume(), expectedScvfArea))
                    DUNE_THROW(Dune::InvalidStateException, "Unexpected scvf geometry area");
                if (localGridGeometry.dofIndex(localGridGeometry.insideScv(scvf)) != insideDofIndex)
                    DUNE_THROW(Dune::InvalidStateException, "Unexpected inside scv");
                if (localGridGeometry.onBoundary(scvf))
                    if (localGridGeometry.numOutsideScvs(scvf) != 0)
                        DUNE_THROW(Dune::InvalidStateException, "NumOutsideScvs != 0 on the boundary");
                if (!localGridGeometry.onBoundary(scvf))
                {
                    // get the normal vector 20 times to mimmick local assembly
                    for (int i = 0; i < numFaceAccesses; ++i)
                    {
                        if (localGridGeometry.numOutsideScvs(scvf) != 1)
                            DUNE_THROW(Dune::InvalidStateException, "NumOutsideScvs != 1 in the domain");
                        const auto& insideScv = localGridGeometry.insideScv(scvf);
                        const auto& outsideScv = localGridGeometry.outsideScv(scvf, 0);
                        const auto n = localGridGeometry.unitOuterNormal(scvf);
                        auto d = localGridGeometry.center(outsideScv) - localGridGeometry.center(insideScv);
                        d /= d.two_norm();
                        if ((d-n).two_norm() > 1e-7)
                            DUNE_THROW(Dune::InvalidStateException, "Unexpected normal/outsideScv alignment");
                    }
                }
            }
            if (scvfCount != numExpectedScvfsPerScv)
                DUNE_THROW(Dune::InvalidStateException, "Unexpected number of scvfs per scv");
        }
    });
}

template<typename Grid>
void testOldGridGeometry(const Grid& grid, std::size_t numCellsPerDirection)
{
    const auto& gridView = grid.leafGridView();
    const double expectedScvfArea = 1.0/numCellsPerDirection;
    const double expectedScvVolume = expectedScvfArea*expectedScvfArea;
    const std::size_t numExpectedScvfsPerScv = 4;
    std::size_t expectedNumScvf = 0;
    for (const auto& e : elements(gridView))
        expectedNumScvf += e.subEntities(1);

    std::cout << "Constructing grid geometry... ";
    auto gridGeometry = printTime("GridGeometry construction", [&] () {
        return Dumux::CCTpfaFVGridGeometry<std::decay_t<decltype(gridView)>, false>{gridView};
    });

    printTime("Running the tests", [&] () {
        if (gridGeometry.numDofs() != gridView.size(0))
            DUNE_THROW(Dune::InvalidStateException, "Unexpected number of dofs");
        if (gridGeometry.numScv() != gridView.size(0))
            DUNE_THROW(Dune::InvalidStateException,
                    "Unexpected number of scvs " << gridGeometry.numScv() << " (expected " << gridView.size(0) << ")");
        if (gridGeometry.numScvf() != expectedNumScvf)
            DUNE_THROW(Dune::InvalidStateException,
                    "Unexpected number of scvfs " << gridGeometry.numScv() << " (expected " << expectedNumScvf << ")");

        auto fvGeometry = localView(gridGeometry);
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            fvGeometry.bind(element);
            for (const auto& scv : scvs(fvGeometry))
            {
                if (!equal(scv.volume(), expectedScvVolume))
                    DUNE_THROW(Dune::InvalidStateException, "Unexpected scv volume");
                if (!equal(scv.geometry().volume(), expectedScvVolume))
                    DUNE_THROW(Dune::InvalidStateException, "Unexpected scv volume");
            }

            int scvfCount = 0;
            for (const auto& scvf : scvfs(fvGeometry))
            {
                scvfCount++;
                if (!equal(scvf.area(), expectedScvfArea))
                    DUNE_THROW(Dune::InvalidStateException, "Unexpected scvf area");
                if (!equal(scvf.geometry().volume(), expectedScvfArea))
                    DUNE_THROW(Dune::InvalidStateException, "Unexpected scvf geometry area");
                if (fvGeometry.scv(scvf.insideScvIdx()).dofIndex() != (*scvs(fvGeometry).begin()).dofIndex())
                    DUNE_THROW(Dune::InvalidStateException, "Unexpected inside scv");
                if (scvf.boundary())
                    if (scvf.numOutsideScvs() > 4)
                        DUNE_THROW(Dune::InvalidStateException, "NumOutsideScvs > 4 on the boundary");
                if (!scvf.boundary())
                {
                    for (int i = 0; i < numFaceAccesses; ++i)
                    {
                        if (scvf.numOutsideScvs() != 1)
                            DUNE_THROW(Dune::InvalidStateException, "NumOutsideScvs != 1 in the domain");
                        const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx(0));
                        const auto n = scvf.unitOuterNormal();
                        auto d = outsideScv.center() - fvGeometry.scv(scvf.insideScvIdx()).center();
                        d /= d.two_norm();
                        if ((d-n).two_norm() > 1e-7)
                            DUNE_THROW(Dune::InvalidStateException, "Unexpected normal/outsideScv alignment");
                    }
                }
            }
            if (scvfCount != numExpectedScvfsPerScv)
                DUNE_THROW(Dune::InvalidStateException, "Unexpected number of scvfs per scv");
        }
    });
}

int main (int argc, char *argv[])
{
    const unsigned int numCells = 2000;
    Dune::YaspGrid<2> grid{
        {1.0, 1.0},
        {numCells, numCells}
    };

    printTime("First run new grid geometry", [&] () {
        testStructuredGrid2D(grid, numCells);
    });

    using namespace std::chrono_literals;
    std::cout << "\nsleeping for 5 seconds so that you can track the memory peaks in you system monitor!\n" << std::endl;
    std::this_thread::sleep_for(5s);

    printTime("First run old grid geometry", [&] () {
        testOldGridGeometry(grid, numCells);
    });

    std::cout << "\nsleeping for 5 seconds so that you can track the memory peaks in you system monitor!\n" << std::endl;
    std::this_thread::sleep_for(5s);

    printTime("second run new grid geometry", [&] () {
        testStructuredGrid2D(grid, numCells);
    });

    std::cout << "\nsleeping for 5 seconds so that you can track the memory peaks in you system monitor!\n" << std::endl;
    std::this_thread::sleep_for(5s);

    printTime("second run old grid geometry", [&] () {
        testOldGridGeometry(grid, numCells);
    });

    std::cout << "\nAll tests passed" << std::endl;
    return 0;
}
