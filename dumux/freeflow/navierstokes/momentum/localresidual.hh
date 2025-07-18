// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesMomentumResidual
 */
#ifndef DUMUX_NAVIERSTOKES_MOMENTUM_LOCAL_RESIDUAL_HH
#define DUMUX_NAVIERSTOKES_MOMENTUM_LOCAL_RESIDUAL_HH

#include <dumux/freeflow/navierstokes/momentum/fcstaggered/localresidual.hh>
#warning "This header is deprecated and will be removed after 3.11. Use NavierStokesMomentumFCStaggeredResidual from dumux/freeflow/navierstokes/momentum/fcstaggered/localresidual.hh instead."

namespace Dumux {
    template<class TypeTag>
    using NavierStokesMomentumResidual = NavierStokesMomentumFCStaggeredResidual<TypeTag>;
}

#endif
