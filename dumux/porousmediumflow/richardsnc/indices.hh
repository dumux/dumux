// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup RichardsNCModel
 * \brief Defines the primary variable and equation indices used by
 *        the richardsnc model.
 */

#ifndef DUMUX_RICHARDSNC_INDICES_HH
#define DUMUX_RICHARDSNC_INDICES_HH

namespace Dumux {

/*!
 * \ingroup RichardsNCModel
 * \brief The indices for the isothermal Richards, n-component model.
 */
struct RichardsNCIndices
{
    //! Component indices
    static constexpr int compMainIdx = 0; //!< main component index

    //! primary variable indices
    static constexpr int pressureIdx = 0; //!< pressure

    //! \note These indices make sense if the first balance is replaced by the
    //!       total mass balance.

    //! Equation indices
    static constexpr int conti0EqIdx = 0; //!< continuity equation index
};

} // end namespace Dumux

#endif
