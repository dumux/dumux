//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/geometry/intersectingentities.hh>

template<int dimworld>
void testIntersectingEntityCartesianGrid()
{
    static constexpr int dim = dimworld;
    using Yasp = Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double, dimworld>>;
    std::array<int, dim> cells; cells.fill(30);
    Dune::FieldVector<double, dimworld> lowerLeft(1.1), upperRight(2.2);
    Yasp yasp(lowerLeft, upperRight, cells);

    for (const auto& element : elements(yasp.leafGridView()))
    {
        const auto center = element.geometry().center();
        // test intersectsPointBoundingBox
        if (!Dumux::intersectsPointBoundingBox(center, lowerLeft, upperRight))
            DUNE_THROW(Dune::Exception, "Element center not in grid bounding box!");

        // test intersectingEntityCartesianGrid
        const auto isIndex = Dumux::intersectingEntityCartesianGrid(center, lowerLeft, upperRight, cells);
        const auto gridIndex = yasp.leafGridView().indexSet().index(element);
        if (isIndex != gridIndex)
            DUNE_THROW(Dune::Exception, "Wrong element index: intersection: " << isIndex
                                        << ", grid index set: " << gridIndex);
    }
}

int main(int argc, char* argv[])
{
    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);

    testIntersectingEntityCartesianGrid<1>();
    testIntersectingEntityCartesianGrid<2>();
    testIntersectingEntityCartesianGrid<3>();

    return 0;
}
