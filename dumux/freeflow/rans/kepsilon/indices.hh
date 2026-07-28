// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \copydoc Dumux::KEpsilonMassIndices
 */
#ifndef DUMUX_RANS_KEPSILON_MASS_INDICES_HH
#define DUMUX_RANS_KEPSILON_MASS_INDICES_HH

#include <dumux/freeflow/navierstokes/mass/1p/indices.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowModels
 * \brief The indices for the two-equation k-epsilon turbulence closure, fused as two additional
 *        transported equations onto the single-phase Navier-Stokes mass balance (pressure
 *        remains equation/index 0, as in NavierStokesMassOnePIndices; k is index 1, epsilon is
 *        index 2).
 */
struct KEpsilonMassIndices : public NavierStokesMassOnePIndices
{
    static constexpr int turbulentKineticEnergyEqIdx = 1;
    static constexpr int turbulentKineticEnergyIdx = turbulentKineticEnergyEqIdx;
    static constexpr int dissipationEqIdx = 2;
    static constexpr int dissipationIdx = dissipationEqIdx;
};

} // end namespace Dumux

#endif
