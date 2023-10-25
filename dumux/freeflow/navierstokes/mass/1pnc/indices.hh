// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesIndices
 */
#ifndef DUMUX_NAVIERSTOKES_MASS_1P_INDICES_HH
#define DUMUX_NAVIERSTOKES_MASS_1P_INDICES_HH

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief The common indices for the isothermal Navier-Stokes mass conservation model.
 */
struct NavierStokesMassOnePIndices
{
    static constexpr int conti0EqIdx = 0; //!< Index of the first (total for pure-fluid systems) mass balance equation
    static constexpr int pressureIdx = conti0EqIdx; //!< Index of the pressure
};

} // end namespace Dumux

#endif
