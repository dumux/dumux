// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterModels
 * \copydoc Dumux::ShallowWaterIndices
 */
#ifndef DUMUX_FREEFLOW_SHALLOW_WATER_INDICES_HH
#define DUMUX_FREEFLOW_SHALLOW_WATER_INDICES_HH

namespace Dumux {

// \{
/*!
 * \ingroup ShallowWaterModels
 * \brief The common indices for the shallow water equations model.
 */
struct ShallowWaterIndices
{
    static constexpr int dimXIdx = 0; //!< Index of the x-component of a vector of size dim
    static constexpr int dimYIdx = 1; //!< Index of the y-component of a vector of size dim

    static constexpr int massBalanceIdx = 0; //!< Index of the mass balance equation
    static constexpr int momentumXBalanceIdx = 1; //!< Index of the x momentum balance equation
    static constexpr int momentumYBalanceIdx = 2; //!< Index of the y momentum balance equation

    static constexpr int waterdepthIdx = massBalanceIdx; //!< Index of the velocity in a solution vector
    static constexpr int velocityXIdx = momentumXBalanceIdx; //!< Index of the velocity in a solution vector
    static constexpr int velocityYIdx = momentumYBalanceIdx; //!< Index of the velocity in a solution vector
    static constexpr int velocityOffset = velocityXIdx; //!< Offset for the velocity index
};

// \}
} // end namespace Dumux

#endif
