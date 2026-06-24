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
 * Model: OnePNC (one phase, two components) with Boussinesq approximation:
 *  - ρ = 1 everywhere in storage and continuity (BoussinesqFluid system)
 *  - D = 1/Ra (Ra appears as inverse diffusion in transport equation)
 *  - Buoyancy density = 1 + C used only in Darcy's law (BoussinesqCVFEDarcyLaw)
 *  - g = (0, -1) dimensionless unit gravity (BoussinesqSpatialParams)
 */
#ifndef DUMUX_BOUSSINESQ_PROPERTIES_HH
#define DUMUX_BOUSSINESQ_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/1pnc/model.hh>

#include "1p_boussinesq_fluidsystem.hh"

#include "boussinesqdarcyslaw.hh"
#include "problem.hh"
#include "spatialparams.hh"

namespace Dumux::Properties {

namespace TTag {
//! Base type tag (discretisation-independent)
struct BoussinesqOneSidedRB { using InheritsFrom = std::tuple<OnePNC>; };
//! Box (CVFE) variant — the only one supported with the Boussinesq Darcy law
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

// Override the Darcy law for Box to use Boussinesq buoyancy (1 + C instead of rho)
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::BoussinesqOneSidedRBBox>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using type = BoussinesqCVFEDarcyLaw<Scalar, GridGeometry>;
};

} // end namespace Dumux::Properties

#endif