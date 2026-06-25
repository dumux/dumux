// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Property definitions for the dimensionless Boussinesq dissolution test (2D and 3D).
 *
 * 2D (BoussinesqOneSidedRB):
 *   - Grid: YaspGrid<2>, gravity in -y
 *   - Primary variables: ψ (streamfunction), C
 *   - numEq = 1 + 1 = 2
 *
 * 3D (BoussinesqOneSidedRB3D):
 *   - Grid: YaspGrid<3>, gravity in -z
 *   - Primary variables: A_x, A_y, A_z (vector potential), C
 *   - numEq = 3 + 1 = 4
 */
#ifndef DUMUX_BOUSSINESQ_PROPERTIES_HH
#define DUMUX_BOUSSINESQ_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include "model/model.hh"
#include "model/modeltraits.hh"
#include "problem.hh"
#include "spatialparams.hh"
#include "1p_boussinesq_fluidsystem.hh"

namespace Dumux::Properties {

// -------------------------------------------------------------------------
// 2D test (existing)
// -------------------------------------------------------------------------
namespace TTag {
struct BoussinesqOneSidedRB
{ using InheritsFrom = std::tuple<BoussinesqVorticityModelBox>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::BoussinesqOneSidedRB>
{ using type = Dune::YaspGrid<2>; };

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::BoussinesqOneSidedRB>
{ using type = BoussinesqVorticityModelTraits<1, 2>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::BoussinesqOneSidedRB>
{ using type = BoussinesqOneSidedRBProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::BoussinesqOneSidedRB>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar       = GetPropType<TypeTag, Properties::Scalar>;
    using type = BoussinesqSpatialParams<GridGeometry, Scalar>;
};

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::BoussinesqOneSidedRB>
{ using type = FluidSystems::BoussinesqFluid<GetPropType<TypeTag, Properties::Scalar>>; };

// -------------------------------------------------------------------------
// 3D test — same dimensionless physics, gravity in -z, vector potential A
// numEq = 3 (A_x, A_y, A_z) + 1 (C) = 4
// -------------------------------------------------------------------------
namespace TTag {
struct BoussinesqOneSidedRB3D
{ using InheritsFrom = std::tuple<BoussinesqVorticityModelBox>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::BoussinesqOneSidedRB3D>
{ using type = Dune::YaspGrid<3>; };

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::BoussinesqOneSidedRB3D>
{ using type = BoussinesqVorticityModelTraits<1, 3>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::BoussinesqOneSidedRB3D>
{ using type = BoussinesqOneSidedRBProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::BoussinesqOneSidedRB3D>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar       = GetPropType<TypeTag, Properties::Scalar>;
    using type = BoussinesqSpatialParams<GridGeometry, Scalar>;
};

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::BoussinesqOneSidedRB3D>
{ using type = FluidSystems::BoussinesqFluid<GetPropType<TypeTag, Properties::Scalar>>; };

} // end namespace Dumux::Properties

#endif
