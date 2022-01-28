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

extern "C" {
#include <stdint.h>
#include <ptscotch.h>
} // end extern C

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

int main (int argc , char **argv)
{
    // if terminal output is desired, e.g. for debugging,
    // set to true, else false
    bool printOutput = true;

    // use MPI helper to initialize MPI
    //
    // This will be needed in order to set up the parallelization and dgraph structure.
    // SCOTCH's dgraphs work analogously to the regular graph, but need some more inputs.
    //
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // TODO add if rank == 0 around everything (read, graph, ...)
    // but consider that load balance is outside of this if,
    // consequently the input params need to be declared outside, too

    // read mesh file
    constexpr int dim = 2; // 3 for ball.msh example
    using Grid = Dune::UGGrid<dim>;
    auto grid = std::shared_ptr<Grid>(Dune::GmshReader<Grid>::read("rectangle.msh"));
    const auto& gridView = grid->leafGridView();
    writeOutput("before_loadbalance", gridView, mpiHelper.rank());

    // create partition array (target rank for each node in graph)
    std::vector<SCOTCH_Num> scotchPartitions;
    if (mpiHelper.rank() == 0)
    {
        using Graph = std::vector<std::vector<int>>;

        // create graph's nodes and edges from mesh
        // (cell-centered methods)
        constexpr int codimension = 0;
        Graph graph(gridView.size(codimension));

        // Displaying element indices
        if (printOutput)
            std::cout << "Element index:" << std::endl;

        for (const auto& element : elements(gridView))
        {
            const auto eIdx = gridView.indexSet().index(element);
            if (printOutput)
                std::cout << eIdx << " ";

            for (const auto& intersection : intersections(gridView, element))
                if (intersection.neighbor())
                    graph[eIdx].push_back( gridView.indexSet().index(intersection.outside()));
        }

        if (printOutput)
            std::cout << std::endl;

        // Number of local graph vertices (cells)
        const SCOTCH_Num numNodes = graph.size();

        // Data structures for graph input to SCOTCH (add 1 for case that
        // graph size is zero)
        std::vector<SCOTCH_Num> vertTab;
        vertTab.reserve(numNodes + 1);
        std::vector<SCOTCH_Num> edgeTab;
        edgeTab.reserve(20*numNodes);

        // Build local graph input for SCOTCH
        // (graph vertices (cells) and
        // number of edges connecting two vertices)
        SCOTCH_Num numEdges = 0;
        vertTab.push_back(0);
        for (auto vertex = graph.begin(); vertex != graph.end(); ++vertex)
        {
            numEdges += vertex->size();
            vertTab.push_back(vertTab.back() + vertex->size());
            edgeTab.insert(edgeTab.end(), vertex->begin(), vertex->end());
        }

        // Shrink vectors to hopefully recover any unused memory
        vertTab.shrink_to_fit();
        edgeTab.shrink_to_fit();

        // initialize graph
        SCOTCH_Graph scotchGraph;
        if (SCOTCH_graphInit(&scotchGraph) != 0)
            DUNE_THROW(Dune::Exception, "Error initializing SCOTCH graph!");

        // check graph's consistency (recommended before building)
        if (SCOTCH_graphCheck(&scotchGraph) != 0)
            DUNE_THROW(Dune::Exception, "Error within SCOTCH graph's consistency!");

        // build graph
        const SCOTCH_Num baseValue = 0;
        if (SCOTCH_graphBuild(&scotchGraph, baseValue, numNodes, vertTab.data(), vertTab.data()+1, NULL, NULL, numEdges, edgeTab.data(), NULL))
            DUNE_THROW(Dune::Exception, "Error building SCOTCH graph!");

        // initialize strategy
        SCOTCH_Strat scotchStrat;
        if (SCOTCH_stratInit(&scotchStrat) != 0)
            DUNE_THROW(Dune::Exception, "Error initializing SCOTCH strategy!");

        // build strategy structure
        // there is also a flag for a "focus" on load balancing called SCOTCH_STRATBALANCE
        // numbers from scotch tests in repo (just to test, if this entire setup works)
        const SCOTCH_Num flagValue = SCOTCH_STRATDEFAULT;
        const SCOTCH_Num numProcessors = mpiHelper.size();
        const double imbalanceRatio = 0.0;

        if (SCOTCH_stratGraphMapBuild(&scotchStrat, flagValue, numProcessors, imbalanceRatio) != 0)
            DUNE_THROW(Dune::Exception, "Error building SCOTCH strategy!");

        // architecture and graphMap are created and called within graphPart
        // if specific ones are desired, one has to call them separately and delete the graphPart function call
        // compute partitioning
        scotchPartitions.resize(numNodes);
        if (SCOTCH_graphPart(&scotchGraph, numProcessors, &scotchStrat, scotchPartitions.data()) != 0)
            DUNE_THROW(Dune::Exception, "Error computing SCOTCH graph mapping!");

        // free memory
        SCOTCH_graphExit(&scotchGraph);
        SCOTCH_stratExit(&scotchStrat);
    }

    // Displaying computed graph partition array
    if (printOutput && !scotchPartitions.empty())
    {
        std::cout << "Graph's partition array:" << std::endl;
        for (int i = 0; i < scotchPartitions.size(); ++i)
            std::cout << scotchPartitions[i] << " ";

        std::cout << std::endl;
    }

    // convert number types
    std::vector<Grid::Rank> partitions(
        scotchPartitions.begin(), scotchPartitions.end()
    );

    // use UG interface to send element to target ranks
    grid->loadBalance(partitions, /*grid level=*/0);

    // write grid out with ranks
    writeOutput("after_loadbalance", gridView, mpiHelper.rank());

    return 0;
}
