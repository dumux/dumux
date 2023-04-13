// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPNCModel
 * \brief Defines the indices required for the two-phase n-component model.
 */

#ifndef DUMUX_2PNC_INDICES_HH
#define DUMUX_2PNC_INDICES_HH

namespace Dumux {

/*!
 * \ingroup TwoPNCModel
 * \brief The indices for the isothermal two-phase n-component model.
 */
struct TwoPNCIndices
{
    // present phases (-> 'pseudo' primary variable)
    static constexpr int firstPhaseOnly = 1;  //!< Only the first phase (in fluid system) is present
    static constexpr int secondPhaseOnly = 2; //!< Only the second phase (in fluid system) is present
    static constexpr int bothPhases = 3;      //!< Both phases are present

    // Primary variable indices
    static constexpr int pressureIdx = 0; //! index for first/second phase pressure (depending on formulation) in privar vector
    static constexpr int switchIdx = 1;   //! index of either the saturation or the mass/mole fraction of the first/second component

    // equation indices
    static constexpr int conti0EqIdx = 0; //! index of the conservation equation for the first component
};

} // end namespace Dumux

#endif
