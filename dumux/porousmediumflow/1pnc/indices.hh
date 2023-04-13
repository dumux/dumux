// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePNCModel
 * \brief Defines the primary variable and equation indices used by
 *        the 1pnc model
 */

#ifndef DUMUX_1PNC_INDICES_HH
#define DUMUX_1PNC_INDICES_HH

namespace Dumux {

/*!
 * \ingroup OnePNCModel
 * \brief The indices for the isothermal one-phase n-component model.
 *
 * \tparam phaseIdx The index of the fluid phase in the fluid system
 */
struct OnePNCIndices
{
    //! Reference index for mass conservation equation.
    static constexpr int conti0EqIdx = 0;
    //! Index for wetting/nonwetting phase pressure (depending on formulation) in a solution vector
    static constexpr int pressureIdx = 0;
};

} // end namespace Dumux

#endif
