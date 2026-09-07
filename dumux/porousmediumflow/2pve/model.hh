// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPVE
 * \brief Adaptation of the fully implicit scheme to the two-phase VE flow model.
 *
 * This model implements two-phase flow of two immiscible fluids
 * \f$\alpha \in \{ w, n \}\f$ using a standard multi-phase Darcy
 * approach as the equation for the conservation of momentum, i.e.
 \f[
 V_\alpha = - \frac{K_{r\alpha}}{\mu_\alpha} \textbf{K}
 \left(\textbf{grad}_{//}\, P_\alpha - \varrho_{\alpha} {\textbf g} \right)
 \f]
 *
 * By inserting this into the equation for the conservation of the
 * phase mass, one gets
 \f[
 \Phi \frac{\partial \varrho_\alpha S_\alpha}{\partial t}
 -
 \text{div} \left\{
 \varrho_\alpha \frac{K_{r\alpha}}{\mu_\alpha} \mathbf{K} \left(\textbf{grad}_{//}\, P_\alpha - \varrho_{\alpha} \mathbf{g} \right)
 \right\} - Q_\alpha = 0 \;,
 \f]
 *
 * All equations are discretized using a cell-centered finite volume scheme as spatial
 * and the implicit Euler method as time discretization.
 *
 * We can reduce the number of unknowns to two by using reconstruction rules for the VE scheme to obtain the non-wetting phase pressure \f$P_n\f$, and relative permeability \f$K_{r\alpha}\f$ and by taking
 * advantage of the fact that \f$S_w + S_n = 1\f$. The capillary pressure can be computed via \f$P_c =
 * P_n - P_w\f$ but this constraint is not required as a closing condition. Currently, the model only supports
 * choosing \f$p_w\f$ and \f$S_n\f$ as primary variables. See \cite Buntic2025 for more details on the model.
 *
 * The current implementation has the following restrictions:
 * - only two- and three-dimensional, structured, axis-aligned grids are
 *   supported
 * - the fine grid must be uniform in the vertical direction
 * - the coarse and fine grids must have identical domain bounds and matching
 *   subdivisions in all horizontal directions
 * - each coarse column must contain exactly one coarse cell in the vertical
 *   direction and the same number of fine cells
 * - each fine-cell center must be contained unambiguously in exactly one
 *   coarse cell
 * - gravity must be nonzero and aligned with the vertical coordinate axis
 * - the wetting phase must be denser than the nonwetting phase
 * - the reconstruction uses a Brooks-Corey material law and does not support
 *   \f$\lambda = 1\f$
 * - only the \f$p_w-S_n\f$ primary-variable formulation is supported
 * - the current upscaling implementation assumes scalar, isotropic
 *   permeability (for computation of coarse-level mobilities, we divide by the permeability)
 *
 * These are restrictions of the current implementation and are not general
 * restrictions of vertical-equilibrium models.
 */

#ifndef DUMUX_TWOPVE_MODEL_HH
#define DUMUX_TWOPVE_MODEL_HH

#include <tuple>

#include <dumux/porousmediumflow/2p/model.hh>
#include <dumux/porousmediumflow/2pve/volumevariables.hh>

namespace Dumux::Properties {

// inherit the complete DuMuX 2p model
namespace TTag {

struct TwoPVE
{
    using InheritsFrom = std::tuple<TwoP>;
};

} // namespace TTag

// only replace the standard 2p volume variables with the VE implementation
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::TwoPVE>
{
private:
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;
    using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;
    using SolidState = GetPropType<TypeTag, Properties::SolidState>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using PermeabilityType = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;
    using DiscretizationMethod = typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod;
    static constexpr bool enableBoxInterfaceSolver = getPropValue<TypeTag, Properties::EnableBoxInterfaceSolver>();
    using SaturationReconstruction = TwoPScvSaturationReconstruction<DiscretizationMethod, enableBoxInterfaceSolver>;

    using Traits = TwoPVolumeVariablesTraits<
        PrimaryVariables,
        FluidSystem,
        FluidState,
        SolidSystem,
        SolidState,
        PermeabilityType,
        ModelTraits,
        SaturationReconstruction>;

public:
    using type = TwoPVEVolumeVariables<Traits>;
};

} // namespace Dumux::Properties

#endif
