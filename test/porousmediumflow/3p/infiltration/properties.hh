// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ThreePTests
 * \brief The properties of the Isothermal NAPL infiltration problem.
 */
#ifndef DUMUX_INFILTRATION_THREEP_PROPERTIES_HH
#define DUMUX_INFILTRATION_THREEP_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/method.hh>
#include <dumux/porousmediumflow/3p/model.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidsystems/1pgas.hh>
#include <dumux/material/fluidsystems/3pimmiscible.hh>
#include <dumux/material/components/tabulatedcomponent.hh>
#include <dumux/material/components/air.hh>
#include <dumux/material/components/mesitylene.hh>
#include <dumux/material/components/h2o.hh>

#include "spatialparams.hh"

#include "problem.hh"

namespace Dumux::Properties {
// Create new type tags
namespace TTag {
struct InfiltrationThreeP { using InheritsFrom = std::tuple<ThreeP>; };
struct InfiltrationThreePBox { using InheritsFrom = std::tuple<InfiltrationThreeP, BoxModel>; };
struct InfiltrationThreePCCTpfa { using InheritsFrom = std::tuple<InfiltrationThreeP, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::InfiltrationThreeP> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::InfiltrationThreeP> { using type = InfiltrationThreePProblem<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::InfiltrationThreeP>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Water = Components::TabulatedComponent<Components::H2O<Scalar>>;
    using WettingFluid = FluidSystems::OnePLiquid<Scalar, Water>;
    using NonwettingFluid = FluidSystems::OnePLiquid<Scalar, Components::Mesitylene<Scalar>>;
    using Gas = FluidSystems::OnePGas<Scalar, Components::Air<Scalar>>;
public:
    using type = FluidSystems::ThreePImmiscible<Scalar, WettingFluid, NonwettingFluid, Gas>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::InfiltrationThreeP>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = InfiltrationThreePSpatialParams<GridGeometry, Scalar>;
};

}// end namespace Dumux::Properties

#endif
