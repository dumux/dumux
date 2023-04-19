// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPOneCModel
 * \copydoc Dumux::TwoPOneCIndices
 */

#ifndef DUMUX_2P1C_INDICES_HH
#define DUMUX_2P1C_INDICES_HH

namespace Dumux {

/*!
 * \ingroup TwoPOneCModel
 * \brief The indices for the two-phase one-component model.
 */
class TwoPOneCIndices
{
public:
    // Present phases (-> 'pseudo' primary variable)
    static const int twoPhases = 1; //!< Both liquid and gas phase are present.
    static const int liquidPhaseOnly = 2; //!< Only the liquid phase is present.
    static const int gasPhaseOnly = 3; //!< Only gas phase is present.

    // Primary variable indices
    static const int pressureIdx = 0; //!< Index for phase pressure in a solution vector.
    static const int switchIdx = 1; //!< Index of saturation or temperature.

    // Equation indices
    static const int conti0EqIdx = 0; //!< Index of the mass conservation equation for the water component.
    static const int energyEqIdx = 1; //<! The index for energy in equation vectors.

    static const int temperatureIdx = -99; //!< For compatibility reasons. Do not use.
};

} // end namespace Dumux

#endif
