// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Nonlinear
 * \brief Updating cached grid variables on a subset of the grid
 */
#ifndef DUMUX_NONLINEAR_PRECONDITIONING_LOCAL_UPDATE_HH
#define DUMUX_NONLINEAR_PRECONDITIONING_LOCAL_UPDATE_HH

#include <type_traits>

#include <dumux/discretization/elementsolution.hh>

namespace Dumux::NonlinearPreconditioning {

/*!
 * \ingroup Nonlinear
 * \brief Bring the cached grid variables up to date on a subset of the elements only
 *
 * A local nonlinear solve changes the solution on one subdomain and must see the effect of that change
 * in the cached volume variables and flux caches. Calling the grid variables' own update would rebuild
 * the caches over the whole grid, which costs more than the local solve it serves and would dominate
 * the run as soon as there are many subdomains.
 *
 * The range must cover the subdomain and its fringe: the fringe elements are never perturbed, but their
 * volume variables are read by the flux terms of the elements on the subdomain boundary.
 *
 * With caching disabled the grid variables hold no state to refresh — volume variables are then built
 * on the fly during assembly — so this reduces to nothing.
 */
template<class GridGeometry, class GridVariables, class SolutionVector, class ElementIndices>
void updateSubdomainVariables(const GridGeometry& gridGeometry,
                              GridVariables& gridVariables,
                              const SolutionVector& sol,
                              const ElementIndices& elementIndices)
{
    auto& gridVolVars = gridVariables.curGridVolVars();
    auto& gridFluxVarsCache = gridVariables.gridFluxVarsCache();

    using GridVolVars = std::decay_t<decltype(gridVolVars)>;
    using GridFluxCache = std::decay_t<decltype(gridFluxVarsCache)>;

    if constexpr (!GridVolVars::cachingEnabled && !GridFluxCache::cachingEnabled)
        return;
    else
    {
        auto fvGeometry = localView(gridGeometry);
        auto elemVolVars = localView(gridVolVars);

        for (const auto eIdx : elementIndices)
        {
            const auto element = gridGeometry.element(eIdx);
            fvGeometry.bindElement(element);

            if constexpr (GridVolVars::cachingEnabled)
            {
                const auto elemSol = elementSolution(element, sol, gridGeometry);
                for (const auto& scv : scvs(fvGeometry))
                    gridVolVars.volVars(scv.dofIndex()).update(elemSol, gridVolVars.problem(), element, scv);
            }
        }

        if constexpr (GridFluxCache::cachingEnabled)
        {
            for (const auto eIdx : elementIndices)
            {
                const auto element = gridGeometry.element(eIdx);
                fvGeometry.bind(element);
                elemVolVars.bind(element, fvGeometry, sol);
                gridFluxVarsCache.updateElement(element, fvGeometry, elemVolVars);
            }
        }
    }
}

} // end namespace Dumux::NonlinearPreconditioning

#endif
