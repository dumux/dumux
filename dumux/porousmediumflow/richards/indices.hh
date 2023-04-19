// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup RichardsModel
 * \brief Index names for the Richards model.
 */

#ifndef DUMUX_RICHARDS_INDICES_HH
#define DUMUX_RICHARDS_INDICES_HH

namespace Dumux {

/*!
 * \ingroup RichardsModel
 * \brief Index names for the Richards model.
 */

struct RichardsIndices
{
    //! Primary variable index for the wetting phase pressure
    static constexpr int pressureIdx = 0;

    //! Equation index for the mass conservation of the wetting phase
    static constexpr int conti0EqIdx = 0;
};

} // end namespace Dumux

#endif
