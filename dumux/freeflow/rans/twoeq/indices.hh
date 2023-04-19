// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoEqModel
 * \copydoc Dumux::RANSTwoEqIndices
 */
#ifndef DUMUX_RANS_TWO_EQ_INDICES_HH
#define DUMUX_RANS_TWO_EQ_INDICES_HH

#include <dumux/freeflow/navierstokes/indices.hh>

namespace Dumux {

// \{
/*!
 * \ingroup TwoEqModel
 * \brief The common indices for isothermal two-equation RANS models.
 *
 * \tparam dimension The dimension of the problem
 * \tparam numComponents The number of considered transported components
 */
template<int dimension, int numComponents>
struct RANSTwoEqIndices : public NavierStokesIndices<dimension>
{
public:
    static constexpr auto turbulentKineticEnergyEqIdx = dimension + numComponents;
    static constexpr auto turbulentKineticEnergyIdx = turbulentKineticEnergyEqIdx;
    static constexpr auto dissipationEqIdx = turbulentKineticEnergyEqIdx + 1;
    static constexpr auto dissipationIdx = dissipationEqIdx;
};

// \}
} // end namespace

#endif
