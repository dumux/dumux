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
#ifndef DUMUX_SOLVERCATEGORY_HH
#define DUMUX_SOLVERCATEGORY_HH

#include <dune/istl/solvers.hh>

namespace Dumux::Detail {

template<class LinearSolverTraits, class GridView>
Dune::SolverCategory::Category solverCategory(const GridView& gridView)
{
    if constexpr (LinearSolverTraits::canCommunicate)
    {
        if (gridView.comm().size() <= 1)
            return Dune::SolverCategory::sequential;

        if (LinearSolverTraits::isNonOverlapping(gridView))
            return Dune::SolverCategory::nonoverlapping;
        else
            return Dune::SolverCategory::overlapping;
    }
    else
        return Dune::SolverCategory::sequential;
}

} // end namespace Dumux::Detail

#endif
