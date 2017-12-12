// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \brief A class for h-adaptivity when using finite volume schemes.
 */
#ifndef DUMUX_ADAPTIVE_MARKELEMENTS_HH
#define DUMUX_ADAPTIVE_MARKELEMENTS_HH

namespace Dumux
{

template<class Grid, class Indicator>
bool markElements(Grid& grid, Indicator& indicator, bool verbose = true)
{
    //! mark elements according to indicator
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

    //! abort if nothing in grid is marked
    const std::size_t sumRefine = grid.comm().sum(refine);
    const std::size_t sumCoarsen = grid.comm().sum(coarsen);

    if (grid.comm().rank() == 0 && verbose)
        std::cout << sumRefine << " cells have been marked to be refined, "
                  << sumCoarsen << " to be coarsened." << std::endl;

    //! return whether or not anything has been marked
    return sumRefine != 0 || sumCoarsen != 0;
}

} //end namespace Dumux

#endif /* DUMUX_ADAPTIVE_MARKELEMENTS_HH */
