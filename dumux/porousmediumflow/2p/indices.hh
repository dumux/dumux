// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPModel
 * \brief Defines the indices required for the two-phase fully implicit model.
 */

#ifndef DUMUX_2P_INDICES_HH
#define DUMUX_2P_INDICES_HH

#include "formulation.hh"

namespace Dumux {

/*!
 * \ingroup TwoPModel
 * \brief Defines the indices required for the two-phase fully implicit model.
 */
struct TwoPIndices
{
    // Primary variable indices
    static constexpr int pressureIdx = 0; //!< index for first/second phase pressure (depending on formulation) in a solution vector
    static constexpr int saturationIdx = 1; //!< index of the saturation of the first/second phase

    // indices of the equations
    static constexpr int conti0EqIdx = 0; //!< index of the first continuity equation
};

} // namespace Dumux


#endif
