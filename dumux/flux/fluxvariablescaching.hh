// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Flux
 * \brief Classes related to flux variables caching
 */
#ifndef DUMUX_DISCRETIZATION_FLUXVAR_CACHING_HH
#define DUMUX_DISCRETIZATION_FLUXVAR_CACHING_HH

namespace Dumux {
namespace FluxVariablesCaching {

//! The empty filler class corresponding to EmptyCache
struct EmptyCacheFiller
{
    EmptyCacheFiller() = default;

    template<class Problem>
    EmptyCacheFiller(const Problem& p) {}

    static constexpr bool isSolDependent = false; // the cache is empty

    template<typename... Args>
    static void fill(Args&&... args) {}
};

//! An empty flux variables cache
template<class S>
struct EmptyCache
{
    //! export type used for scalar values
    using Scalar = S;

    template<typename... Args>
    void update(Args&&... args) {}
};

#ifndef DOXYGEN // hide the empty caches from doxygen
// an empty cache filler
// \note Never use the _EmptyCache directly as it lead to ambiguous definitions
struct _EmptyCache
{ using Filler = EmptyCacheFiller; };
#endif // DOXYGEN

/*!
 * \ingroup Flux
 * \brief Empty caches to use in a constitutive flux law/process, e.g. Darcy's law
 */
struct EmptyAdvectionCache : public _EmptyCache {};
struct EmptyDiffusionCache : public _EmptyCache {};
struct EmptyHeatConductionCache : public _EmptyCache {};

} // end namespace FluxVariablesCaching
} // end namespace Dumux

#endif
