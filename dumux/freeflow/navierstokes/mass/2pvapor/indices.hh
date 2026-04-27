// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_FREEFLOW_NAVIERSTOKES_MASS_2PVAPOR_INDICES_HH
#define DUMUX_FREEFLOW_NAVIERSTOKES_MASS_2PVAPOR_INDICES_HH

#include <dumux/freeflow/navierstokes/mass/2p/indices.hh>

namespace Dumux {

/*!
 * \brief Extends the 2p Cahn-Hilliard indices with a vapor transport equation.
 *
 * Primary variables: [p, φ, μ, c_v]
 *   c_v [kg/m³]: vapor mass concentration in the gas phase
 */
struct NavierStokesMassTwoPVaporIndices : public NavierStokesMassTwoPIndices
{
    static constexpr int vaporEqIdx = 3; //!< Index of the vapor transport equation
    static constexpr int vaporIdx   = vaporEqIdx; //!< Index of the vapor concentration
};

} // end namespace Dumux

#endif
