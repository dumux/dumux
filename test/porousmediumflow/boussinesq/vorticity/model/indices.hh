// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_BOUSSINESQ_VORTICITY_INDICES_HH
#define DUMUX_BOUSSINESQ_VORTICITY_INDICES_HH

namespace Dumux {

/*!
 * \brief Indices for the vector-potential Boussinesq model.
 *
 * Primary variables layout:
 *   [0 .. numPotentialEqs-1]          vector potential components A_k
 *   [numPotentialEqs .. numEq-1]      scalar transport components C_i
 *
 * In 2D only the z-component of the vector potential is non-trivial
 * (the classical streamfunction ψ), so numPotentialEqs = 1.
 * In 3D all dim components are solved, numPotentialEqs = dim.
 *
 * \tparam nComp  number of transported scalar components (≥ 1)
 * \tparam dim    spatial dimension (2 or 3)
 */
template<int nComp, int dim>
struct BoussinesqVorticityIndices
{
    static constexpr int numPotentialEqs = (dim == 2) ? 1 : dim;

    // ---- vector-potential (A_k) ----
    static constexpr int vectorPotentialIdx(int k)    { return k; }
    static constexpr int vectorPotentialEqIdx(int k)  { return k; }

    // ---- scalar transport (C_i) ----
    static constexpr int transportIdx(int i)          { return numPotentialEqs + i; }
    static constexpr int transportEqIdx(int i)        { return numPotentialEqs + i; }

    // ---- backward-compat aliases for single-component 2D usage ----
    static constexpr int streamfunctionIdx   = 0;
    static constexpr int streamfunctionEqIdx = 0;
    static constexpr int concentrationIdx    = numPotentialEqs;   // == transportIdx(0)
    static constexpr int concentrationEqIdx  = numPotentialEqs;
};

} // namespace Dumux

#endif
