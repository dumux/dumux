// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 *
 * \copydoc Dumux::NavierStokesVolumeVariables
 */
#ifndef DUMUX_NAVIERSTOKES_MOMENTUM_PQ1BUBBLE_VOLUME_VARIABLES_HH
#define DUMUX_NAVIERSTOKES_MOMENTUM_PQ1BUBBLE_VOLUME_VARIABLES_HH

#warning "This header is deprecated and will be removed after 3.7"

#include <dumux/freeflow/navierstokes/momentum/cvfe/volumevariables.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Volume variables for the single-phase Navier-Stokes model.
 */
template <class Traits>
using NavierStokesMomentumPQ1BubbleVolumeVariables [[deprecated("Use NavierStokesMomentumCVFEVolumeVariables")]] = NavierStokesMomentumCVFEVolumeVariables<Traits>;

} // end namespace Dumux

#endif
