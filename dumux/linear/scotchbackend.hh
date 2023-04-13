// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

/*!
 * \file
 * \ingroup Linear
 * \brief An interface to the scotch library for matrix reordering
 * \note You need to have PTSCOTCH installed to use this feature
 */
#ifndef DUMUX_SCOTCH_BACKEND_HH
#define DUMUX_SCOTCH_BACKEND_HH

#include <string>
#include <vector>
#include <iostream>

#include <dune/common/exceptions.hh>

#if HAVE_PTSCOTCH
#if HAVE_MPI
#include <mpi.h>
#endif
extern "C"
{
#include <stdint.h>
#include <ptscotch.h>
}
#else
#warning "PTSCOTCH was not found on your system. Dumux::ScotchBackend won't do anything."
#endif

namespace Dumux {

#if HAVE_PTSCOTCH

/*!
 * \ingroup Linear
 * \brief A wrapper around a SCOTCH graph object
 */
template<class IndexType = int>
class ScotchGraph
{
public:
    //! a graph represented by an adjacency list
    using Graph = std::vector<std::vector<IndexType>>;

    ScotchGraph(const Graph& graph)
    {
        // Number of local graph vertices
        const SCOTCH_Num numNodes = graph.size();

        // Data structures for graph input to SCOTCH
        // add 1 for case that graph size is zero
        vertTab_.reserve(numNodes + 1);
        edgeTab_.reserve(20*numNodes);

        // Build local graph input for SCOTCH
        // (number of graph vertices (cells),
        //  number of edges connecting two vertices)
        SCOTCH_Num numEdges = 0;
        vertTab_.push_back(0);
        for (auto vertex = graph.begin(); vertex != graph.end(); ++vertex)
        {
            numEdges += vertex->size();
            vertTab_.push_back(vertTab_.back() + vertex->size());
            edgeTab_.insert(edgeTab_.end(), vertex->begin(), vertex->end());
        }

        // Shrink vectors to hopefully recover any unused memory
        vertTab_.shrink_to_fit();
        edgeTab_.shrink_to_fit();

        if (SCOTCH_graphInit(&scotchGraph_) != 0)
            DUNE_THROW(Dune::Exception, "Error initializing SCOTCH graph!");

        // check graph's consistency (recommended before building)
        if (SCOTCH_graphCheck(&scotchGraph_) != 0)
            DUNE_THROW(Dune::Exception, "Error within SCOTCH graph's consistency!");

        // build graph
        const SCOTCH_Num baseValue = 0; // C-style array indexing
        if (SCOTCH_graphBuild(&scotchGraph_, baseValue, numNodes, vertTab_.data(), vertTab_.data()+1, NULL, NULL, numEdges, edgeTab_.data(), NULL))
            DUNE_THROW(Dune::Exception, "Error building SCOTCH graph!");
    }

    //! Clean-up the graph
    ~ScotchGraph()
    {
        SCOTCH_graphExit(&scotchGraph_);
    }

    //! Get the raw point to the data (to pass to C interface)
    SCOTCH_Graph* data()
    { return &scotchGraph_; }

private:
    SCOTCH_Graph scotchGraph_;
    // we have to maintain these ourselves to keep the Scotch graph valid
    std::vector<SCOTCH_Num> vertTab_;
    std::vector<SCOTCH_Num> edgeTab_;
};

/*!
 * \ingroup Linear
 * \brief A wrapper around a SCOTCH strategy object
 */
class ScotchGraphOrderStrategy
{
public:
    ScotchGraphOrderStrategy(const std::string& strategy = "")
    {
        if (SCOTCH_stratInit(&strategy_) != 0)
            DUNE_THROW(Dune::Exception, "Error initializing SCOTCH strategy!");

        // Set SCOTCH strategy (if provided)
        if (!strategy.empty())
            SCOTCH_stratGraphOrder(&strategy_, strategy.c_str());
    }

    //! Clean-up the strategy
    ~ScotchGraphOrderStrategy()
    {
        SCOTCH_stratExit(&strategy_);
    }

    //! Get the raw point to the data (to pass to C interface)
    SCOTCH_Strat* data()
    { return &strategy_; }

private:
    SCOTCH_Strat strategy_;
};

#endif // HAVE_PTSCOTCH

/*!
 * \ingroup Linear
 * \brief A reordering backend using the scotch library
 * \note You need to have PTSCOTCH installed to use this feature
 */
template<class IndexType>
class ScotchBackend
{
public:
    //! the graph type
    using Graph = std::vector<std::vector<IndexType>>;

    //! Compute reordering (map[old] -> new) using
    //! Gibbs-Poole-Stockmeyer (GPS) re-ordering
    static std::vector<int> computeGPSReordering(const Graph& graph,
                                                 std::size_t numPasses = 5)
    {
        // Create strategy string for Gibbs-Poole-Stockmeyer ordering
        std::string strategy = "g{pass= " + std::to_string(numPasses) + "}";
        return computeReordering(graph, strategy);
    }

    //! Compute graph re-ordering
    static std::vector<int> computeReordering(const Graph& graph,
                                              std::string scotchStrategy = "")
    {
        std::vector<int> permutation, inversePermutation;
        computeReordering(graph, permutation, inversePermutation, scotchStrategy);
        return permutation;
    }

    //! Compute graph re-ordering
    static void computeReordering(const Graph& graph,
                                   std::vector<int>& permutation,
                                   std::vector<int>& inversePermutation,
                                   std::string scotchStrategy = "")
    {
#if HAVE_PTSCOTCH
        ScotchGraph<IndexType> scotchGraph(graph);
        ScotchGraphOrderStrategy strategy(scotchStrategy);

        // Vector to hold permutation vectors
        const auto graphSize = graph.size();
        std::vector<SCOTCH_Num> permutationIndices(graphSize);
        std::vector<SCOTCH_Num> inversePermutationIndices(graphSize);

        // Reset SCOTCH random number generator to produce deterministic partitions
        SCOTCH_randomReset();

        // Compute re-ordering
        if (SCOTCH_graphOrder(
            scotchGraph.data(), strategy.data(), permutationIndices.data(),
            inversePermutationIndices.data(), NULL, NULL, NULL
        ))
            DUNE_THROW(Dune::Exception, "Error reordering SCOTCH graph!");

        // Copy permutation vectors
        permutation.resize(graphSize);
        inversePermutation.resize(graphSize);
        std::copy(permutationIndices.begin(), permutationIndices.end(),
                  permutation.begin());
        std::copy(inversePermutationIndices.begin(), inversePermutationIndices.end(),
                  inversePermutation.begin());
#endif // HAVE_PTSCOTCH
    }
};

} // end namespace Dumux

#endif
