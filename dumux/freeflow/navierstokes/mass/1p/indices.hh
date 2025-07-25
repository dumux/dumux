// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesMassOnePIndices
 */
#ifndef DUMUX_FREEFLOW_NAVIERSTOKES_MASS_1P_INDICES_HH
#define DUMUX_FREEFLOW_NAVIERSTOKES_MASS_1P_INDICES_HH

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief The common indices for the isothermal Navier-Stokes mass conservation model.
 */
struct NavierStokesMassOnePIndices
{
    static constexpr int conti0EqIdx = 0; //!< Index of the mass/mole balance equation
    static constexpr int pressureIdx = conti0EqIdx; //!< Index of the pressure
};

} // end namespace Dumux

#endif
