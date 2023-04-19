// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FluidSystems
 * \brief @copybrief Dumux::NullParameterCache
 */
#ifndef DUMUX_NULL_PARAMETER_CACHE_HH
#define DUMUX_NULL_PARAMETER_CACHE_HH

#include "parametercachebase.hh"

namespace Dumux {
/*!
 * \ingroup FluidSystems
 * \brief The a parameter cache which does nothing
 */
class NullParameterCache : public ParameterCacheBase<NullParameterCache> {};

} // end namespace

#endif
