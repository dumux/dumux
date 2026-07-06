// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Property definitions for the pressure-based dimensionless Boussinesq
 *        dissolution test (2D, one-sided Rayleigh-Bénard).
 *
 * Model: DuMux OnePNC with Box discretisation.
 *   - 2 components (solvent idx 0, solute idx 1)
 *   - replaceCompEqIdx = 0  →  eq 0 = total mass balance (∇·u = 0, drives p)
 *   - eq 1 = solute transport
 *   - AdvectionType = BoussinesqCVFEDarcyLaw  →  ρ_buoy = 1 + C in gravity term
 *   - BoussinesqFluid  →  ρ = 1 everywhere (Boussinesq storage/transport)
 *   - useMoles = false  →  primary var 1 = mass fraction C ∈ [0,1]
 *
 * Pressure reference fixed by Robin BC at top:
 *   neumann[0] = (p - p_ref) * penalty
 */
#ifndef DUMUX_BOUSSINESQ_PRESSURE_PROPERTIES_HH
#define DUMUX_BOUSSINESQ_PRESSURE_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/1pnc/model.hh>

#include "problem.hh"
#include "../spatialparams.hh"
#include "../1p_boussinesq_fluidsystem.hh"
#include "../boussinesqdarcyslaw.hh"

namespace Dumux::Properties {

namespace TTag {
struct BoussinesqPressureTest
{ using InheritsFrom = std::tuple<OnePNC, BoxModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::BoussinesqPressureTest>
{ using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<double, 2>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::BoussinesqPressureTest>
{ using type = BoussinesqPressureProblem<TypeTag>; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::BoussinesqPressureTest>
{ using type = FluidSystems::BoussinesqFluid<GetPropType<TypeTag, Properties::Scalar>>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::BoussinesqPressureTest>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar       = GetPropType<TypeTag, Properties::Scalar>;
    using type = BoussinesqSpatialParams<GridGeometry, Scalar>;
};

// Replace default CVFE Darcy law with Boussinesq version (buoyancy = 1 + C)
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::BoussinesqPressureTest>
{
    using Scalar       = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using type = BoussinesqCVFEDarcyLaw<Scalar, GridGeometry>;
};

// Use mass fractions: primary var 1 = C ∈ [0,1] directly
template<class TypeTag>
struct UseMoles<TypeTag, TTag::BoussinesqPressureTest>
{ static constexpr bool value = false; };

// Replace the solvent (comp 0) equation with the total mass balance.
// → eq 0 becomes ∇·u = 0, the elliptic pressure equation.
template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::BoussinesqPressureTest>
{ static constexpr int value = 0; };

} // end namespace Dumux::Properties

#endif
