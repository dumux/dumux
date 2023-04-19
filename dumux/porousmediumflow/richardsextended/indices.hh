// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ExtendedRichardsModel
 * \brief Index names for the extended Richards model.
 */

#ifndef DUMUX_RICHARDSEXTENDED_INDICES_HH
#define DUMUX_RICHARDSEXTENDED_INDICES_HH

#include <dumux/porousmediumflow/richards/indices.hh>

namespace Dumux {

/*!
 * \ingroup ExtendedRichardsModel
 * \brief Index names for the extended Richards model.
 */

struct ExtendedRichardsIndices : public RichardsIndices
{
    static constexpr int switchIdx = 0;

    // present phases (-> 'pseudo' primary variable)
    static constexpr int liquidPhaseOnly = 1; //!< Only the liquid phase is present
    static constexpr int gasPhaseOnly = 2; //!< Only the gas phase is present
    static constexpr int bothPhases = 3; //!< Both phases are present
};

} // end namespace Dumux

#endif
