// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CahnHilliardModel
 * \brief Defines the indices for the Cahn-Hilliard model
 */
#ifndef DUMUX_PHASEFIELD_CAHNHILLIARD_INDICES_HH
#define DUMUX_PHASEFIELD_CAHNHILLIARD_INDICES_HH

namespace Dumux {

/*!
 * \ingroup CahnHilliardModel
 * \brief The indices for the Cahn-Hilliard model.
 */
struct CahnHilliardIndices
{
    // primary variable indices
    static constexpr int concentrationIdx = 0;
    static constexpr int chemicalPotentialIdx = 1;

    // equation indices
    static constexpr int massBalanceEqIdx = 0;
    static constexpr int chemicalPotentialEqIdx = 1;
};

} // end namespace Dumux

#endif
