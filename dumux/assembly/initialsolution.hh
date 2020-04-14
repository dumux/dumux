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
 * \ingroup Assembly
 * \brief Function to create initial solution vectors
 */
#ifndef DUMUX_ASSEMBLY_INITIAL_SOLUTION_HH
#define DUMUX_ASSEMBLY_INITIAL_SOLUTION_HH

#include <vector>
#include <type_traits>

#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \file
 * \ingroup Assembly
 * \brief Set a solution vector to the initial solution provided by the problem
 */
template<class SolutionVector, class Problem>
void assembleInitialSolution(SolutionVector& sol, const Problem& problem)
{
    const auto& gg = problem.gridGeometry();
    using GridGeometry = std::decay_t<decltype(gg)>;

    // box method
    if constexpr (GridGeometry::discMethod == DiscretizationMethod::box)
    {
        constexpr int dim = GridGeometry::GridView::dimension;
        const auto numDofs = gg.vertexMapper().size();
        const auto numVert = gg.gridView().size(dim);
        sol.resize(numDofs);

        // if there are more dofs than vertices (enriched nodal dofs), we have to
        // call initial for all dofs at the nodes, coming from all neighboring elements.
        if (numDofs != numVert)
        {
            std::vector<bool> dofVisited(numDofs, false);
            for (const auto& element : elements(gg.gridView()))
            {
                for (int i = 0; i < element.subEntities(dim); ++i)
                {
                    const auto dofIdxGlobal = gg.vertexMapper().subIndex(element, i, dim);
                    // forward to implementation if value at dof is not set yet
                    if (!dofVisited[dofIdxGlobal])
                    {
                        sol[dofIdxGlobal] = problem.initial(element.template subEntity<dim>(i));
                        dofVisited[dofIdxGlobal] = true;
                    }
                }
            }
        }

        // otherwise we directly loop over the vertices
        else
        {
            for (const auto& vertex : vertices(gg.gridView()))
                sol[gg.vertexMapper().index(vertex)] = problem.initial(vertex);
        }
    }

    // staggered methods
    else if constexpr (GridGeometry::discMethod == DiscretizationMethod::staggered)
    {
        problem.applyInitialSolution(sol);
    }

    // default: cell-centered methods
    else
    {
        sol.resize(gg.numDofs());
        for (const auto& element : elements(gg.gridView()))
            sol[gg.elementMapper().index(element)] = problem.initial(element);
    }
}

/*!
 * \file
 * \ingroup Assembly
 * \brief Create a solution vector filled with the initial solution provided by the problem
 */
template<class SolutionVector, class Problem>
SolutionVector makeInitialSolution(const Problem& problem)
{
    SolutionVector sol;
    assembleInitialSolution(sol, problem);
    return sol;
}

} // end namespace Dumux

#endif
