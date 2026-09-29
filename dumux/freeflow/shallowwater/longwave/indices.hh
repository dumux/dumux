// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup LongWaveModel
 * \copydoc Dumux::LongWaveIndices
 */
#ifndef DUMUX_FREEFLOW_SHALLOWWATER_LONGWAVE_INDICES_HH
#define DUMUX_FREEFLOW_SHALLOWWATER_LONGWAVE_INDICES_HH

namespace Dumux {

/*!
 * \ingroup LongWaveModel
 * \brief The indices for the long-wave approximations of the shallow water equations
 */
struct LongWaveIndices
{
    static constexpr int massBalanceIdx = 0; //!< Index of the mass balance equation
    static constexpr int waterDepthIdx = 0; //!< Index of the water depth in a solution vector
};

} // end namespace Dumux

#endif
