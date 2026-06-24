// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Property definitions for the dimensionless Boussinesq dissolution test.
 *
 * Model: OnePNC (one phase, two components) with streamfunction-Boussinesq formulation:
 *  - Primary variable 0: ψ (streamfunction, in the pressure slot)
 *  - Primary variable 1: C (solute mass fraction)
 *  - Eq 0: ∇·(∇ψ − C·êₓ) = 0  (streamfunction Poisson, no time derivative)
 *  - Eq 1: φ ∂C/∂t + ∇·(u C − D∇C) = 0  with u = (∂ψ/∂y, −∂ψ/∂x)
 *  - D = 1/Ra (from BoussinesqFluid system)
 *  - BCs: ψ=0 Dirichlet on all walls; C=1 Dirichlet at top; ∂C/∂n=0 elsewhere
 */
#ifndef DUMUX_BOUSSINESQ_PROPERTIES_HH
#define DUMUX_BOUSSINESQ_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/1pnc/model.hh>

#include "1p_boussinesq_fluidsystem.hh"
#include "localresidual.hh"
#include "streamfunctionadvection.hh"
#include "problem.hh"
#include "spatialparams.hh"

namespace Dumux::Properties {

namespace TTag {
struct BoussinesqOneSidedRB { using InheritsFrom = std::tuple<OnePNC>; };
struct BoussinesqOneSidedRBBox
{ using InheritsFrom = std::tuple<BoussinesqOneSidedRB, BoxModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::BoussinesqOneSidedRB>
{ using type = Dune::YaspGrid<2>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::BoussinesqOneSidedRB>
{ using type = BoussinesqOneSidedRBProblem<TypeTag>; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::BoussinesqOneSidedRB>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::BoussinesqFluid<Scalar>;
};

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::BoussinesqOneSidedRB>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = BoussinesqSpatialParams<GridGeometry, Scalar>;
};

// Use mass fractions. With M_solvent = M_solute = 1 this is equivalent to mole fractions.
template<class TypeTag>
struct UseMoles<TypeTag, TTag::BoussinesqOneSidedRB>
{ static constexpr bool value = false; };

// Equation 0 is the streamfunction Poisson equation (repurposed from total mass balance).
// Setting replaceCompEqIdx = 0 gives the right model-traits structure (2 eqs, eq 0 without
// a component-storage term in the parent, which we then zero out in our local residual).
template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::BoussinesqOneSidedRB>
{ static constexpr int value = 0; };

// Streamfunction-Boussinesq local residual: overrides eq 0 storage (zero) and flux.
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::BoussinesqOneSidedRB>
{ using type = BoussinesqStreamfunctionLocalResidual<TypeTag>; };

// Advection from curl(ψ) instead of Darcy's law.
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::BoussinesqOneSidedRBBox>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using type = StreamfunctionAdvection<Scalar, GridGeometry>;
};

} // end namespace Dumux::Properties

#endif
