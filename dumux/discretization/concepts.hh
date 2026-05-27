// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Concepts for discretization types
 */
#ifndef DUMUX_DISCRETIZATION_CONCEPTS_HH
#define DUMUX_DISCRETIZATION_CONCEPTS_HH

namespace Dumux::Experimental::Concepts {

/*!
 * \ingroup Discretization
 * \brief Concept for finite-volume discretizations.
 *
 * A discretization satisfies this concept if it provides iteration
 * over sub-control volumes and sub-control volume faces via the functions
 * `scvs()` and `scvfs()`.
 */
template<class ED>
concept FVElementDiscretization = requires(const ED& ed) {
    scvs(ed);
    scvfs(ed);
};

/*!
 * \ingroup Discretization
 * \brief Concept for hybrid finite-element/finite-volume discretizations.
 *
 * Satisfies `FVElementDiscretization` and additionally provides `nonCVLocalDofs()`
 * for the degrees of freedom which are not related to control volumes.
 */
template<class ED>
concept HybridElementDiscretization = FVElementDiscretization<ED>
    && requires(const ED& ed) { nonCVLocalDofs(ed); };

/*!
 * \ingroup Discretization
 * \brief Concept for pure finite-element discretizations (no FV structure).
 *
 * Provides `nonCVLocalDofs()` but does NOT satisfy
 * `FVElementDiscretization`.
 */
template<class ED>
concept FEElementDiscretization = (!FVElementDiscretization<ED>)
    && requires(const ED& ed) { nonCVLocalDofs(ed); };

} // namespace Dumux::Experimental::Concepts

#endif // DUMUX_DISCRETIZATION_CONCEPTS_HH
