// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesMassTwoPIndices
 */
#ifndef DUMUX_FREEFLOW_NAVIERSTOKES_MASS_2P_INDICES_HH
#define DUMUX_FREEFLOW_NAVIERSTOKES_MASS_2P_INDICES_HH

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief The common indices for the isothermal Navier-Stokes mass conservation model.
 */
struct NavierStokesMassTwoPIndices
{
    static constexpr int conti0EqIdx = 0; //!< Index of the mass balance equation
    static constexpr int pressureIdx = conti0EqIdx; //!< Index of the pressure

    // two more equations for the Cahn-Hilliard type phase field model
    // the phase field transport equation and the chemical potential equation

    static constexpr int phaseFieldEqIdx = 1; //!< Index phase field transport balance equation
    static constexpr int phaseFieldIdx = phaseFieldEqIdx; //!< Index of the phase field

    static constexpr int chemicalPotentialEqIdx = 2; //!< Index of the chemical potential equation
    static constexpr int chemicalPotentialIdx = chemicalPotentialEqIdx; //!< Index of the chemical potential
};

} // end namespace Dumux

#endif
