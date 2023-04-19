// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ThreePModel
 * \brief Defines the indices for the three-phase model.
 */

#ifndef DUMUX_3P_INDICES_HH
#define DUMUX_3P_INDICES_HH

namespace Dumux {

/*!
 * \ingroup ThreePModel
 * \brief The common indices for the isothermal three-phase model.
 */
class ThreePIndices
{
public:
    // Primary variable indices
    static constexpr int pressureIdx = 0; //!< index for gas phase pressure in a solution vector
    static constexpr int swIdx = 1; //!< index of water (more wetting than the other liquid) saturation
    static constexpr int snIdx = 2; //!< index of (e.g.) NAPL saturation

    // equation indices
    static constexpr int conti0EqIdx = 0; //!< index of first balance equation
};

} // end namespace Dumux

#endif
