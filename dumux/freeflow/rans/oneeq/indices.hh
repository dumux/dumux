// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OneEqModel
 * \copydoc Dumux::OneEqIndices
 */
#ifndef DUMUX_ONEEQ_INDICES_HH
#define DUMUX_ONEEQ_INDICES_HH

#include <dumux/freeflow/navierstokes/indices.hh>

namespace Dumux {

// \{
/*!
 * \ingroup OneEqModel
 * \brief The common indices for the isothermal one-equation turbulence model by Spalart-Allmaras
 *
 * \tparam dimension The dimension of the problem
 * \tparam numComponents The number of considered transported components
 */
template<int dimension, int numComponents>
struct OneEqIndices : public NavierStokesIndices<dimension>
{
public:
    static constexpr auto viscosityTildeEqIdx = dimension + numComponents;
    static constexpr auto viscosityTildeIdx = viscosityTildeEqIdx;
};

// \}
} // end namespace

#endif
