// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup HeleShawModel
 * \brief Defines the indices for the hybrid CVFE Hele-Shaw two-phase model
 */
#ifndef DUMUX_PHASEFIELD_HELESHAW_2P_CVFE_INDICES_HH
#define DUMUX_PHASEFIELD_HELESHAW_2P_CVFE_INDICES_HH

namespace Dumux {

/*!
 * \ingroup HeleShawModel
 * \brief The indices for the hybrid CVFE Hele-Shaw two-phase model.
 */
struct HeleShawTwoPCVFEIndices
{
    // primary variable indices
    static constexpr int pressureIdx   = 0;
    static constexpr int phaseFieldIdx = 1;
    static constexpr int chemPotIdx    = 2;

    // equation indices (1:1 with primary variables here)
    static constexpr int continuityEqIdx = 0;
    static constexpr int phaseFieldEqIdx = 1;
    static constexpr int chemPotEqIdx    = 2;
};

} // end namespace Dumux

#endif
