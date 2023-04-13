// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ThreePWaterOilTests
 * \brief The properties of the non-isothermal steam-assisted gravity
 *        drainage (SAGD) problem.
 */

#ifndef DUMUX_SAGDPROPERTIES_HH
#define DUMUX_SAGDPROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/3pwateroil/model.hh>

#include <dumux/material/fluidsystems/h2oheavyoil.hh>
#include <dumux/material/solidsystems/1csolid.hh>
#include <dumux/material/components/constant.hh>

#include "problem.hh"
#include "spatialparams.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct Sagd { using InheritsFrom = std::tuple<ThreePWaterOilNI>; };
struct ThreePWaterOilSagdBox { using InheritsFrom = std::tuple<Sagd, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::Sagd> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Sagd> { using type = Dumux::SagdProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Sagd>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = SagdSpatialParams<GridGeometry, Scalar>;
};

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Sagd>
{ using type = Dumux::FluidSystems::H2OHeavyOil<GetPropType<TypeTag, Properties::Scalar>>; };

template<class TypeTag>
struct OnlyGasPhaseCanDisappear<TypeTag, TTag::Sagd> { static constexpr bool value = true; };

template<class TypeTag>
struct UseMoles<TypeTag, TTag::Sagd> { static constexpr bool value = true; };

// Set the solid system
template<class TypeTag>
struct SolidSystem<TypeTag, TTag::Sagd>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using InertComponent = Components::Constant<1, Scalar>;
    using type = SolidSystems::InertSolidPhase<Scalar, InertComponent>;
};

} // end namespace Dumux::Properties

#endif
