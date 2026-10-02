// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test the element coloring of the partial reassembler for the cell-centered TPFA scheme
 */
#include <config.h>

#include <array>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <dumux/common/initialize.hh>
#include <dumux/assembly/entitycolor.hh>
#include <dumux/assembly/partialreassembler.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>

namespace Dumux {

// the partial reassembler only queries the grid geometry of the assembler for the coloring
template<class GG>
class GridGeometryOnlyAssembler
{
public:
    using Scalar = double;
    using GridGeometry = GG;
    using JacobianMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<Scalar, 1, 1>>;

    GridGeometryOnlyAssembler(const GridGeometry& gridGeometry)
    : gridGeometry_(gridGeometry)
    {}

    const GridGeometry& gridGeometry() const
    { return gridGeometry_; }

private:
    const GridGeometry& gridGeometry_;
};

// elements that changed, and all elements sharing a face with one of them
template<class GridGeometry>
std::set<std::size_t> expectedRedElements(const GridGeometry& gridGeometry,
                                          const std::set<std::size_t>& changed)
{
    std::set<std::size_t> red(changed);
    const auto& mapper = gridGeometry.elementMapper();
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        const auto eIdx = mapper.index(element);
        for (const auto& intersection : intersections(gridGeometry.gridView(), element))
            if (intersection.neighbor() && changed.count(mapper.index(intersection.outside())))
                red.insert(eIdx);
    }
    return red;
}

template<class GridGeometry>
bool checkColors(const GridGeometry& gridGeometry,
                 const std::set<std::size_t>& changed,
                 const std::string& testCase)
{
    using Assembler = GridGeometryOnlyAssembler<GridGeometry>;
    const Assembler assembler(gridGeometry);
    PartialReassembler<Assembler> reassembler(assembler);

    const double threshold = 1e-6;
    const auto numElements = gridGeometry.elementMapper().size();
    std::vector<double> distance(numElements, 0.1*threshold);
    for (const auto eIdx : changed)
        distance[eIdx] = 10.0*threshold;

    reassembler.computeColors(assembler, distance, threshold);

    const auto red = expectedRedElements(gridGeometry, changed);
    std::size_t numWrong = 0;
    for (std::size_t eIdx = 0; eIdx < numElements; ++eIdx)
    {
        const auto expected = red.count(eIdx) ? EntityColor::red : EntityColor::green;
        if (reassembler.elementColor(eIdx) != expected || reassembler.dofColor(eIdx) != expected)
        {
            ++numWrong;
            std::cout << testCase << ": element " << eIdx << " should be "
                      << (expected == EntityColor::red ? "red" : "green") << std::endl;
        }
    }

    if (numWrong > 0)
    {
        std::cout << testCase << ": " << numWrong << " of " << numElements
                  << " elements have the wrong color" << std::endl;
        return false;
    }

    std::cout << testCase << ": " << red.size() << " red elements, as expected" << std::endl;
    return true;
}

} // end namespace Dumux

int main(int argc, char* argv[])
{
    using namespace Dumux;
    initialize(argc, argv);

    using Grid = Dune::YaspGrid<2>;
    using GridGeometry = CCTpfaFVGridGeometry<typename Grid::LeafGridView>;

    constexpr int numCellsPerDirection = 5;
    const Grid grid(Dune::FieldVector<double, 2>(1.0),
                    std::array<int, 2>{{numCellsPerDirection, numCellsPerDirection}});
    const GridGeometry gridGeometry(grid.leafGridView());

    // YaspGrid numbers elements lexicographically
    const auto index = [](int i, int j) -> std::size_t { return j*numCellsPerDirection + i; };

    bool passed = true;
    passed &= checkColors(gridGeometry, {}, "no change");
    passed &= checkColors(gridGeometry, {index(2, 2)}, "interior element");
    passed &= checkColors(gridGeometry, {index(0, 0)}, "corner element");
    passed &= checkColors(gridGeometry, {index(1, 1), index(3, 3)}, "two separate elements");
    passed &= checkColors(gridGeometry, {index(2, 2), index(3, 2)}, "two adjacent elements");

    return passed ? 0 : 1;
}
