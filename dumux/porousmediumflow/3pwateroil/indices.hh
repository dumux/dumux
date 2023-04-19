// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ThreePWaterOilModel
 * \brief Defines the indices required for the 3p2cni model.
 */

#ifndef DUMUX_3P2CNI_INDICES_HH
#define DUMUX_3P2CNI_INDICES_HH

namespace Dumux {

/*!
 * \ingroup ThreePWaterOilModel
 * \brief The indices for the isothermal 3p2cni model.
 */
class ThreePWaterOilIndices
{
public:
    // present phases (-> 'pseudo' primary variable)
    static const int threePhases = 1; //!< All three phases are present
    static const int wPhaseOnly = 2; //!< Only the water phase is present
    static const int gnPhaseOnly = 3; //!< Only gas and NAPL phases are present
    static const int wnPhaseOnly = 4; //!< Only water and NAPL phases are present
    static const int gPhaseOnly = 5; //!< Only gas phase is present
    static const int wgPhaseOnly = 6; //!< Only water and gas phases are present

    // Primary variable indices
    static const int pressureIdx = 0; //!< Index for phase pressure in a solution vector
    static const int switch1Idx = 1; //!< Index 1 of saturation, mole fraction or temperature
    static const int switch2Idx = 2; //!< Index 2 of saturation, mole fraction or temperature

    // equation indices
    static const int conti0EqIdx = 0; //!< Index of the mass conservation equation for the water component
};

} // end namespace Dumux

#endif
