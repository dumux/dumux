// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \copydoc Dumux::OneEqMassIndices
 */
#ifndef DUMUX_RANS_ONEEQ_MASS_INDICES_HH
#define DUMUX_RANS_ONEEQ_MASS_INDICES_HH

#include <dumux/freeflow/navierstokes/mass/1p/indices.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowModels
 * \brief The indices for the one-equation (Spalart-Allmaras) turbulence closure, fused as a
 *        second transported equation onto the single-phase Navier-Stokes mass balance
 *        (pressure remains equation/index 0, exactly as in NavierStokesMassOnePIndices;
 *        the working viscosity ν̃ is equation/index 1).
 */
struct OneEqMassIndices : public NavierStokesMassOnePIndices
{
    static constexpr int viscosityTildeEqIdx = 1;
    static constexpr int viscosityTildeIdx = viscosityTildeEqIdx;
};

} // end namespace Dumux

#endif
