// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test the parallel index set of the linear solvers on a grid with a process without elements
 */
#include <config.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <dune/common/bigunsignedint.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/owneroverlapcopy.hh>

#include <dumux/common/initialize.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/parallelhelpers.hh>
#include <dumux/linear/solvercategory.hh>

namespace Dumux::Test {

//! Every degree of freedom has exactly one owner, and the copies of the other processes receive its value
template<class GridGeometry>
void checkParallelIndexSet(const GridGeometry& gridGeometry, std::size_t numGlobalDofs, const std::string& name)
{
    using Traits = LinearSolverTraits<GridGeometry>;
    using Comm = Dune::OwnerOverlapCopyCommunication<Dune::bigunsignedint<96>, int>;

    const auto& gridView = gridGeometry.gridView();
    Comm comm(gridView.comm(), Detail::solverCategory<Traits>(gridView));
    const ParallelISTLHelper<Traits> helper(gridView, gridGeometry.dofMapper());
    helper.createParallelIndexSet(comm);

    std::size_t numOwned = 0;
    Dune::BlockVector<Dune::FieldVector<double, 1>> ownerRank(gridGeometry.numDofs());
    for (std::size_t i = 0; i < ownerRank.size(); ++i)
    {
        numOwned += helper.isOwned(i);
        ownerRank[i] = helper.isOwned(i) ? gridView.comm().rank() : -1.0;
    }
    comm.copyOwnerToAll(ownerRank, ownerRank);

    numOwned = gridView.comm().sum(numOwned);
    if (numOwned != numGlobalDofs)
        DUNE_THROW(Dune::Exception, name << ": " << numOwned << " owned degrees of freedom, expected " << numGlobalDofs);

    for (std::size_t i = 0; i < ownerRank.size(); ++i)
        if (ownerRank[i][0] < 0.0)
            DUNE_THROW(Dune::Exception, name << ": degree of freedom " << i << " on process " << gridView.comm().rank()
                                        << " has not received the value of its owner");
}

} // end namespace Dumux::Test

int main(int argc, char** argv)
{
    using namespace Dumux;

    initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    constexpr int dim = 2;
    using Grid = Dune::UGGrid<dim>;
    const std::array<unsigned int, dim> cells = {8, 8};
    auto grid = Dune::StructuredGridFactory<Grid>::createCubeGrid({0.0, 0.0}, {1.0, 1.0}, cells);
    const auto gridView = grid->leafGridView();

    // the last process gets no elements, so with more than two processes the others share degrees of freedom
    // and the last one shares none
    std::vector<Grid::Rank> partitions;
    if (mpiHelper.rank() == 0)
    {
        const std::size_t numElements = gridView.size(0);
        const std::size_t numTargets = std::max(mpiHelper.size() - 1, 1);
        partitions.resize(numElements);
        for (std::size_t eIdx = 0; eIdx < numElements; ++eIdx)
            partitions[eIdx] = static_cast<Grid::Rank>(eIdx*numTargets/numElements);
    }
    grid->loadBalance(partitions, /*grid level=*/0);

    std::size_t numInteriorElements = 0;
    for ([[maybe_unused]] const auto& element : elements(gridView, Dune::Partitions::interior))
        ++numInteriorElements;
    if (mpiHelper.size() > 1 && mpiHelper.rank() == mpiHelper.size() - 1 && numInteriorElements > 0)
        DUNE_THROW(Dune::Exception, "The last process holds " << numInteriorElements << " elements");

    Test::checkParallelIndexSet(CCTpfaFVGridGeometry<Grid::LeafGridView>(gridView), cells[0]*cells[1], "cell-centered");
    Test::checkParallelIndexSet(BoxFVGridGeometry<double, Grid::LeafGridView>(gridView), (cells[0] + 1)*(cells[1] + 1), "box");

    if (mpiHelper.rank() == 0)
        std::cout << "Parallel index set is consistent on " << mpiHelper.size() << " processes" << std::endl;

    return 0;
}
