// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Adaptive
 * \brief A function to mark element for refinement or coarsening
 */
#ifndef DUMUX_ADAPTIVE_MARKELEMENTS_HH
#define DUMUX_ADAPTIVE_MARKELEMENTS_HH

#include <dune/grid/common/rangegenerators.hh>

namespace Dumux {

/*!
 * \ingroup Adaptive
 * \brief A function to mark element for refinement or coarsening
 * \param grid the grid to mark the entities on
 * \param indicator indicated per element whether to refine, coarsen, do nothing. It has a ()-operator taking an element
 * \param verbose If verbose output to std::cout is enabled
 * \return bool whether or not anything has been marked
 */
template<class Grid, class Indicator>
bool markElements(Grid& grid, const Indicator& indicator, bool verbose = true)
{
    // mark elements according to indicator
    std::size_t refine = 0;
    std::size_t coarsen = 0;
    for (const auto& element : elements(grid.leafGridView(), Dune::Partitions::interior))
    {
        const auto mark = indicator(element);
        if (mark > 0)
        {
            grid.mark(1,  element);
            ++refine;
        }
        else if (mark < 0)
        {
            grid.mark(-1, element);
            ++coarsen;
        }
    }

    // abort if nothing in grid is marked
    const std::size_t sumRefine = grid.comm().sum(refine);
    const std::size_t sumCoarsen = grid.comm().sum(coarsen);

    if (grid.comm().rank() == 0 && verbose)
        std::cout << sumRefine << " cells have been marked to be refined, "
                  << sumCoarsen << " to be coarsened." << std::endl;

    // return whether or not anything has been marked
    return sumRefine != 0 || sumCoarsen != 0;
}

} //end namespace Dumux

#endif
