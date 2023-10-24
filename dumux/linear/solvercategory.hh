// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Linear
 * \brief Solver category
 */
#ifndef DUMUX_LINEAR_SOLVERCATEGORY_HH
#define DUMUX_LINEAR_SOLVERCATEGORY_HH

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
    {
        if (gridView.comm().size() > 1)
            DUNE_THROW(Dune::InvalidStateException,
                "Attempt to construct parallel solver but LinearSolverTraits::canCommunicate is false. " <<
                "Maybe the grid implementation does not support distributed parallelism."
            );
    }

    return Dune::SolverCategory::sequential;
}

} // end namespace Dumux::Detail

#endif
