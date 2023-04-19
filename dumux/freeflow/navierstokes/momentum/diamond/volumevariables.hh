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
#ifndef DUMUX_NAVIERSTOKES_MOMENTUM_DIAMOND_VOLUME_VARIABLES_HH
#define DUMUX_NAVIERSTOKES_MOMENTUM_DIAMOND_VOLUME_VARIABLES_HH

#warning "This header is deprecated and will be removed after 3.7"

#include <dumux/freeflow/navierstokes/momentum/cvfe/volumevariables.hh>

namespace Dumux {

template<class TypeTag>
using NavierStokesMomentumDiamondVolumeVariables [[deprecated("Use NavierStokesMomentumCVFEVolumeVariables")]] = NavierStokesMomentumCVFEVolumeVariables<TypeTag>;

} // end namespace Dumux

#endif
