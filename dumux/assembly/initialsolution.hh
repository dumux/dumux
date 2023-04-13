// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
 * \ingroup Assembly
 * \brief Set a solution vector to the initial solution provided by the problem
 */
template<class SolutionVector, class Problem>
void assembleInitialSolution(SolutionVector& sol, const Problem& problem)
{
    const auto& gg = problem.gridGeometry();
    using GridGeometry = std::decay_t<decltype(gg)>;

    // box method
    if constexpr (GridGeometry::discMethod == DiscretizationMethods::box)
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
    else if constexpr (GridGeometry::discMethod == DiscretizationMethods::staggered)
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
