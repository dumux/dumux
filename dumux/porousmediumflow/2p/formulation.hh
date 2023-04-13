// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPModel
 * \brief Defines an enumeration for the formulations accepted by the two-phase model.
 */

#ifndef DUMUX_2P_FORMULATION_INDICES_HH
#define DUMUX_2P_FORMULATION_INDICES_HH

namespace Dumux {

/*!
 * \ingroup TwoPModel
 * \brief Enumerates the formulations which the two-phase model accepts.
 */
enum class TwoPFormulation
{
    p0s1, //!< first phase pressure and second phase saturation as primary variables
    p1s0  //!< first phase saturation and second phase pressure as primary variables
};

} // end namespace Dumux

#endif
