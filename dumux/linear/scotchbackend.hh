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
 * \ingroup Linear
 * \brief An interface to the scotch library for matrix reordering
 * \note You need to have PTSCOTCH installed to use this feature
 */
#ifndef DUMUX_SCOTCH_BACKEND_HH
#define DUMUX_SCOTCH_BACKEND_HH

#include <string>
#include <vector>
#include <iostream>

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

/*!
 * \file
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
        // Number of local graph vertices (cells)
        const SCOTCH_Num vertnbr = graph.size();

        // Data structures for graph input to SCOTCH (add 1 for case that
        // graph size is zero)
        std::vector<SCOTCH_Num> verttab;
        verttab.reserve(vertnbr + 1);
        std::vector<SCOTCH_Num> edgetab;
        edgetab.reserve(20*vertnbr);

        // Build local graph input for SCOTCH
        // (number of graph vertices (cells),
        //  number of edges connecting two vertices)
        SCOTCH_Num edgenbr = 0;
        verttab.push_back(0);
        typename Graph::const_iterator vertex;
        for (vertex = graph.begin(); vertex != graph.end(); ++vertex)
        {
            edgenbr += vertex->size();
            verttab.push_back(verttab.back() + vertex->size());
            edgetab.insert(edgetab.end(), vertex->begin(), vertex->end());
        }

        // Shrink vectors to hopefully recover an unused memory
        verttab.shrink_to_fit();
        edgetab.shrink_to_fit();

        // Create SCOTCH graph
        SCOTCH_Graph scotchGraph;

        // C-style array indexing
        const SCOTCH_Num baseval = 0;

        // Create SCOTCH graph and initialise
        if (SCOTCH_graphInit(&scotchGraph) != 0)
        {
            std::cerr << "Error initializing SCOTCH graph!" << std::endl;
            exit(1);
        }

        // Build SCOTCH graph
        if (SCOTCH_graphBuild(&scotchGraph, baseval,
                              vertnbr, &verttab[0], &verttab[1], NULL, NULL,
                              edgenbr, &edgetab[0], NULL))
        {
            std::cerr << "Error building SCOTCH graph!" << std::endl;
            exit(1);
        }

        // Re-ordering strategy
        SCOTCH_Strat strat;
        SCOTCH_stratInit(&strat);

        // Set SCOTCH strategy (if provided)
        if (!scotchStrategy.empty())
        SCOTCH_stratGraphOrder(&strat, scotchStrategy.c_str());

        // Vector to hold permutation vectors
        std::vector<SCOTCH_Num> permutationIndices(vertnbr);
        std::vector<SCOTCH_Num> inversePermutationIndices(vertnbr);

        // Reset SCOTCH random number generator to produce deterministic partitions
        SCOTCH_randomReset();

        // Compute re-ordering
        if (SCOTCH_graphOrder(&scotchGraph, &strat, permutationIndices.data(),
                              inversePermutationIndices.data(), NULL, NULL, NULL))
        {
            std::cerr << "SCOTCH: Error during reordering graph!" << std::endl;
            exit(1);
        }

        // Clean up SCOTCH objects
        SCOTCH_graphExit(&scotchGraph);
        SCOTCH_stratExit(&strat);

        // Copy permutation vectors
        permutation.resize(vertnbr);
        inversePermutation.resize(vertnbr);
        std::copy(permutationIndices.begin(), permutationIndices.end(),
                  permutation.begin());
        std::copy(inversePermutationIndices.begin(), inversePermutationIndices.end(),
                  inversePermutation.begin());
#endif // HAVE_PTSCOTCH
    }
};

} // end namespace Dumux

#endif
