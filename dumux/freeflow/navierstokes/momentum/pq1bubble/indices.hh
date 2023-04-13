// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesIndices
 */
#ifndef DUMUX_NAVIERSTOKES_MOMENTUM_PQ1BUBBLE_INDICES_HH
#define DUMUX_NAVIERSTOKES_MOMENTUM_PQ1BUBBLE_INDICES_HH


#warning "This header is deprecated and will be removed after 3.7"

#include <dumux/freeflow/navierstokes/momentum/cvfe/indices.hh>

namespace Dumux {

template <int dimension>
using NavierStokesMomentumPQ1BubbleIndices [[deprecated("Use NavierStokesMomentumCVFEIndices")]] = NavierStokesMomentumCVFEIndices<dimension>;

} // end namespace Dumux

#endif
