// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_BOUSSINESQ_VORTICITY_MODEL_TRAITS_HH
#define DUMUX_BOUSSINESQ_VORTICITY_MODEL_TRAITS_HH

#include "indices.hh"

namespace Dumux {

/*!
 * \brief Model traits for the vector-potential Boussinesq formulation.
 *
 * Equation layout:
 *   dim == 2:  1  vector-potential (streamfunction ψ) + nComp transport
 *   dim == 3:  3  vector-potential components (A_x, A_y, A_z) + nComp transport
 *
 * \tparam nComp  number of transported scalar components (default 1)
 * \tparam dim    spatial dimension (default 2)
 */
template<int nComp = 1, int dim = 2>
struct BoussinesqVorticityModelTraits
{
    using Indices = BoussinesqVorticityIndices<nComp, dim>;

    static constexpr int numPotentialEqs      = (dim == 2) ? 1 : dim;
    static constexpr int numEq()              { return numPotentialEqs + nComp; }
    static constexpr int numFluidPhases()     { return 1; }
    static constexpr int numFluidComponents() { return nComp; }

    static constexpr bool enableAdvection()          { return true; }
    static constexpr bool enableMolecularDiffusion() { return true; }
    static constexpr bool enableEnergyBalance()      { return false; }
};

} // namespace Dumux

#endif
