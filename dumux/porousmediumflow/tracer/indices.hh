// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TracerModel
 * \brief Defines the primary variable and equation indices used by the isothermal tracer model.
 */

#ifndef DUMUX_TRACER_INDICES_HH
#define DUMUX_TRACER_INDICES_HH

namespace Dumux {

// \{

/*!
 * \ingroup TracerModel
 * \brief Defines the primary variable and equation indices used by the isothermal tracer model.
 */
struct TracerIndices
{
    /*!
     * Component indices are just numbered by component index
     * primary variable indices are just numbered by component index
     * Equation indices
     */
    static const int transportEqIdx = 0; //!< transport equation index
};

// \}
}

#endif
