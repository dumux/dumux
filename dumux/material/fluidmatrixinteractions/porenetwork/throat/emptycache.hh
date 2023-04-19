// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Fluidmatrixinteractions
 * \ingroup PoreNetworkModels
 * \brief An empty cache for transmissibility laws using only standard quantities
 */
#ifndef DUMUX_PNM_EMPTY_CACHE_HH
#define DUMUX_PNM_EMPTY_CACHE_HH

namespace Dumux::PoreNetwork {

struct EmptyCache
{
    template<class ...Args>
    void fill(Args&&...) {}
};

} // end namespace Dumux::PoreNetwork

#endif
