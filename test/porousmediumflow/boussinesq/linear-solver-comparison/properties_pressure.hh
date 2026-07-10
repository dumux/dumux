// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Property definitions for the pressure-based Boussinesq model, used
 *        by the linear-solver comparison benchmark (small grid).
 *
 * Identical TypeTag/physics to white-noise-perturbations/properties_pressure.hh
 * -- duplicated here rather than included so this directory can keep its own
 * small-grid params.input independent of the production tests.
 */
#ifndef DUMUX_BOUSSINESQ_LINSOLVECOMPARE_PRESSURE_PROPERTIES_HH
#define DUMUX_BOUSSINESQ_LINSOLVECOMPARE_PRESSURE_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/1pnc/model.hh>

#include "../pressure/problem.hh"
#include "../pressure/boussinesqdarcyslaw.hh"
#include "../common/spatialparams.hh"
#include "../common/1p_boussinesq_fluidsystem.hh"

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
