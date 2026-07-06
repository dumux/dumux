// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup AllenCahnModel
 * \brief Defines the indices for the Allen-Cahn model
 */
#ifndef DUMUX_PHASEFIELD_ALLENCAHN_INDICES_HH
#define DUMUX_PHASEFIELD_ALLENCAHN_INDICES_HH

namespace Dumux {

/*!
 * \ingroup AllenCahnModel
 * \brief The indices for the Allen-Cahn model.
 */
struct AllenCahnIndices
{
    // primary variable indices
    static constexpr int phaseFieldIdx = 0;

    // equation indices
    static constexpr int phaseFieldEqIdx = 0;
};

} // end namespace Dumux

#endif
