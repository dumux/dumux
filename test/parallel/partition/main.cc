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

#include <config.h>

#include <memory>
#include <vector>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dumux/parallel/scotchpartitioner.hh>

template<class GV>
void writeOutput(const std::string& name, const GV& gridView, int myRank)
{
    Dune::VTKWriter<GV> writer(gridView);
    std::vector<int> rank(gridView.size(0)); // resizes to number of elements (codim-0-entities) on this process only
    for (const auto& element : elements(gridView, Dune::Partitions::interior))
    {
        const auto eIdx = gridView.indexSet().index(element);
        rank[eIdx] = myRank;
    }
    writer.addCellData(rank, "rank");
    writer.write(name);
}

template<class GridView>
typename Dumux::ScotchPartitioner<typename GridView::Grid::Rank>::Graph
buildCCTpfaGraph(const GridView& gridView)
{
    constexpr int codimension = 0;
    using Graph = typename Dumux::ScotchPartitioner<typename GridView::Grid::Rank>::Graph;
    Graph graph(gridView.size(codimension));

    for (const auto& element : elements(gridView))
    {
        const auto eIdx = gridView.indexSet().index(element);
        for (const auto& intersection : intersections(gridView, element))
            if (intersection.neighbor())
                graph[eIdx].push_back(gridView.indexSet().index(intersection.outside()));
    }

    return graph;
}

void test2D(int rank, std::size_t numProcessors)
{
    constexpr int dim = 2;
    using Grid = Dune::UGGrid<dim>;
    auto grid = std::shared_ptr<Grid>(Dune::GmshReader<Grid>::read("rectangle.msh"));
    const auto& gridView = grid->leafGridView();
    writeOutput("before_loadbalance_2d", gridView, rank);

    // create partition array (target rank for each node in graph)
    std::vector<Grid::Rank> partitions;
    if (rank == 0)
    {
        const auto graph = buildCCTpfaGraph(gridView);
        Dumux::ScotchPartitioner<Grid::Rank>::partition(graph, numProcessors, partitions);

        if (numProcessors == 2 || numProcessors == 4)
            if (std::count(partitions.begin(), partitions.end(), 0) != std::count(partitions.begin(), partitions.end(), 1))
                DUNE_THROW(Dune::Exception, "Unexpected unbalanced partition!");
    }

    grid->loadBalance(partitions, /*grid level=*/0);
    writeOutput("after_loadbalance_2d", gridView, rank);
}

void test3D(int rank, std::size_t numProcessors)
{
    constexpr int dim = 3;
    using Grid = Dune::UGGrid<dim>;
    auto grid = std::shared_ptr<Grid>(Dune::GmshReader<Grid>::read("ball.msh"));
    const auto& gridView = grid->leafGridView();
    writeOutput("before_loadbalance_3d", gridView, rank);

    // create partition array (target rank for each node in graph)
    std::vector<Grid::Rank> partitions;
    if (rank == 0)
    {
        const auto graph = buildCCTpfaGraph(gridView);
        Dumux::ScotchPartitioner<Grid::Rank>::partition(graph, numProcessors, partitions);

        if (numProcessors == 2 || numProcessors == 4)
            if (std::count(partitions.begin(), partitions.end(), 0) != std::count(partitions.begin(), partitions.end(), 1))
                DUNE_THROW(Dune::Exception, "Unexpected unbalanced partition!");
    }

    grid->loadBalance(partitions, /*grid level=*/0);
    writeOutput("after_loadbalance_3d", gridView, rank);
}

int main (int argc , char **argv)
{
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    test2D(mpiHelper.rank(), mpiHelper.size());
    test3D(mpiHelper.rank(), mpiHelper.size());

    return 0;
}
