// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Test that every degree of freedom on a periodic boundary is mapped to all its periodic images
 */
#include <config.h>

#include <algorithm>
#include <array>
#include <bitset>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/spgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/deprecated.hh>
#include <dumux/io/grid/gridmanager_sp.hh>

#include <dumux/discretization/box/fvgridgeometry.hh>
#include <dumux/discretization/pq1bubble/fvgridgeometry.hh>
#include <dumux/discretization/pq2/fvgridgeometry.hh>
#include <dumux/discretization/pq3/fvgridgeometry.hh>
#include <dumux/discretization/pq2/fegriddiscretization.hh>
#include <dumux/discretization/pq3/fegriddiscretization.hh>

namespace Dumux {

template<class GridGeometry>
auto dofPositions(const GridGeometry& gridGeometry)
{
    using GlobalPosition = typename GridGeometry::GridView::template Codim<0>::Geometry::GlobalCoordinate;
    std::vector<GlobalPosition> positions(gridGeometry.numDofs());
    auto localGeometry = localView(gridGeometry);
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        localGeometry.bindElement(element);
        for (const auto& localDof : localDofs(localGeometry))
            positions[localDof.dofIndex()] = ipData(localGeometry, localDof).global();
    }
    return positions;
}

template<class GridGeometry, std::size_t dim>
void checkPeriodicDofMap(const GridGeometry& gridGeometry, const std::bitset<dim>& periodic, const std::string& name)
{
    const auto positions = dofPositions(gridGeometry);
    const auto bBoxMin = gridGeometry.bBoxMin();
    const auto bBoxMax = gridGeometry.bBoxMax();
    const double eps = 1e-8;

    std::map<std::size_t, std::size_t> groupSizeCount;
    for (std::size_t dofIdx = 0; dofIdx < gridGeometry.numDofs(); ++dofIdx)
    {
        const auto& pos = positions[dofIdx];
        int numPeriodicDirections = 0;
        for (int dir = 0; dir < dim; ++dir)
            if (periodic[dir] && std::min(pos[dir] - bBoxMin[dir], bBoxMax[dir] - pos[dir]) < eps)
                ++numPeriodicDirections;

        if ((numPeriodicDirections > 0) != gridGeometry.dofOnPeriodicBoundary(dofIdx))
            DUNE_THROW(Dune::Exception, name << ": dof " << dofIdx << " at " << pos
                       << " is wrongly " << (numPeriodicDirections > 0 ? "not " : "") << "on a periodic boundary");

        if (numPeriodicDirections == 0)
            continue;

        std::vector<std::size_t> group{ dofIdx };
        for (const auto periodicDofIdx : Deprecated::rangeOfPeriodicallyMappedDofs(gridGeometry, dofIdx))
        {
            const auto& periodicPos = positions[periodicDofIdx];
            for (int dir = 0; dir < dim; ++dir)
            {
                const auto distance = std::abs(periodicPos[dir] - pos[dir]);
                const bool isImage = distance < eps
                    || (periodic[dir] && std::abs(distance - (bBoxMax[dir] - bBoxMin[dir])) < eps);
                if (!isImage)
                    DUNE_THROW(Dune::Exception, name << ": dof " << dofIdx << " at " << pos
                               << " is mapped to dof " << periodicDofIdx << " at " << periodicPos
                               << ", which is not a periodic image");
            }
            group.push_back(periodicDofIdx);
        }

        std::ranges::sort(group);
        const std::size_t expectedGroupSize = std::size_t(1) << numPeriodicDirections;
        if (std::ranges::adjacent_find(group) != group.end() || group.size() != expectedGroupSize)
            DUNE_THROW(Dune::Exception, name << ": dof " << dofIdx << " at " << pos
                       << " is mapped to " << group.size() - 1 << " distinct dofs instead of " << expectedGroupSize - 1);

        ++groupSizeCount[group.size()];
    }

    std::cout << name << ":";
    for (const auto& [size, count] : groupSizeCount)
        std::cout << " " << count << " dofs in groups of " << size << ";";
    std::cout << std::endl;
}

template<int dim>
void testAllSchemes(const Dune::FieldVector<double, dim>& upperRight,
                    const std::array<int, dim>& cells,
                    const std::bitset<dim>& periodic)
{
    using Grid = Dune::SPGrid<double, dim>;
    GridManager<Grid> gridManager;
    gridManager.init(Dune::FieldVector<double, dim>(0.0), upperRight, cells, "", 1, periodic);
    const auto& gridView = gridManager.grid().leafGridView();
    using GridView = std::decay_t<decltype(gridView)>;

    const std::string suffix = " (" + std::to_string(dim) + "d, periodic " + periodic.to_string() + ")";
    checkPeriodicDofMap(BoxFVGridGeometry<double, GridView, true>(gridView), periodic, "box" + suffix);
    checkPeriodicDofMap(PQ1BubbleFVGridGeometry<double, GridView, true>(gridView), periodic, "pq1bubble" + suffix);
    checkPeriodicDofMap(PQ2FVGridGeometry<double, GridView, true>(gridView), periodic, "pq2" + suffix);
    checkPeriodicDofMap(PQ3FVGridGeometry<double, GridView, true>(gridView), periodic, "pq3" + suffix);
    checkPeriodicDofMap(Experimental::PQ2FEGridDiscretization<double, GridView>(gridView), periodic, "pq2 fe" + suffix);
    checkPeriodicDofMap(Experimental::PQ3FEGridDiscretization<double, GridView>(gridView), periodic, "pq3 fe" + suffix);
}

} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;

    initialize(argc, argv);
    Parameters::init();

    testAllSchemes<2>({2.0, 1.0}, {4, 3}, std::bitset<2>("11"));
    testAllSchemes<2>({2.0, 1.0}, {4, 3}, std::bitset<2>("01"));
    testAllSchemes<3>({1.0, 2.0, 1.5}, {3, 3, 3}, std::bitset<3>("111"));

    return 0;
}
