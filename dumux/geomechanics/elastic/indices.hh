// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Elastic
 * \brief Defines the indices for the elastic model
 */
#ifndef DUMUX_ELASTIC_INDICES_HH
#define DUMUX_ELASTIC_INDICES_HH

namespace Dumux {

/*!
 * \ingroup Elastic
 * \brief The indices for the linear elasticity model.
 */
struct ElasticIndices
{
    // returns the equation index for a given space direction
    static constexpr int momentum(int dirIdx) { return dirIdx; };

    // returns the primary variable index for a given space direction
    static constexpr int u(int dirIdx) { return dirIdx; };

    // Equation indices
    static constexpr int momentumXEqIdx = 0;
    static constexpr int momentumYEqIdx = 1;
    static constexpr int momentumZEqIdx = 2;

    // primary variable indices
    static constexpr int uxIdx = 0;
    static constexpr int uyIdx = 1;
    static constexpr int uzIdx = 2;
};

} // end namespace Dumux

#endif
