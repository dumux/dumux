// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePModel
 * \brief  Defines the indices for the one-phase fully implicit model.
 */

#ifndef DUMUX_1P_INDICES_HH
#define DUMUX_1P_INDICES_HH

namespace Dumux {
// \{

/*!
 * \ingroup OnePModel
 * \brief Indices for the one-phase model.
 *
 * \tparam offset The first index in a primary variable vector.
 */
template<int offset = 0>
struct OnePIndices
{
    static const int PVOffset = offset;      //!< the first index in primary variable vector
    static const int conti0EqIdx = PVOffset; //!< index for the mass balance
    static const int pressureIdx = PVOffset; //!< index of the primary variable
};

// \}
} // end namespace

#endif
