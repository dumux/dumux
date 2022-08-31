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
 * \ingroup Parallel
 * \brief An interface to the Scotch library for graph partitioning
 * \note You need to have PTSCOTCH installed to use this feature
 */
#ifndef DUMUX_PARALLEL_SCOTCH_PARTITIONER_HH
#define DUMUX_PARALLEL_SCOTCH_PARTITIONER_HH

#include <string>
#include <vector>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dumux/linear/scotchbackend.hh>

namespace Dumux {

#if HAVE_PTSCOTCH

/*!
 * \ingroup Parallel
 * \brief A wrapper around a SCOTCH strategy object
 */
class ScotchGraphMapStrategy
{
public:
    ScotchGraphMapStrategy(std::size_t numProcessors, double imbalanceRatio = 0.0, int flag = SCOTCH_STRATDEFAULT)
    {
        if (SCOTCH_stratInit(&strategy_) != 0)
            DUNE_THROW(Dune::Exception, "Error initializing SCOTCH strategy!");

        if (SCOTCH_stratGraphMapBuild(&strategy_,  static_cast<SCOTCH_Num>(flag), static_cast<SCOTCH_Num>(numProcessors), imbalanceRatio) != 0)
            DUNE_THROW(Dune::Exception, "Error building SCOTCH strategy!");
    }

    //! Clean-up the strategy
    ~ScotchGraphMapStrategy()
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
 * \ingroup Parallel
 * \brief A reordering backend using the scotch library
 * \note You need to have PTSCOTCH installed to use this feature
 */
template<class IndexType = int>
class ScotchPartitioner
{
public:
    //! the graph type
    using Graph = std::vector<std::vector<IndexType>>;

    //! Compute graph partition
    static std::vector<IndexType> partition(const Graph& graph, std::size_t numProcessors)
    {
        std::vector<IndexType> targetProcessors;
        partition(graph, numProcessors, targetProcessors);
        return targetProcessors;
    }

    //! Compute graph partition
    static void partition(const Graph& graph, std::size_t numProcessors,
                          std::vector<IndexType>& targetProcessors)
    {
#if HAVE_PTSCOTCH
        ScotchGraph<IndexType> scotchGraph(graph);
        ScotchGraphMapStrategy strategy(numProcessors);

        // architecture and graphMap are created and called within graphPart
        // if specific ones are desired, one has to call them separately and delete the graphPart function call
        // compute partitioning
        const auto graphSize = graph.size();
        std::vector<SCOTCH_Num> scotchPartitions(graphSize);
        if (SCOTCH_graphPart(scotchGraph.data(), static_cast<SCOTCH_Num>(numProcessors), strategy.data(), scotchPartitions.data()) != 0)
            DUNE_THROW(Dune::Exception, "Error computing SCOTCH graph mapping!");

        // convert number types
        targetProcessors = std::vector<IndexType>(
            scotchPartitions.begin(), scotchPartitions.end()
        );
#endif // HAVE_PTSCOTCH
    }
};

} // end namespace Dumux

#endif
