// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#ifndef DUMUX_MICP_COLUMN_SIMPLE_CHEM_PROPERTIES_HH
#define DUMUX_MICP_COLUMN_SIMPLE_CHEM_PROPERTIES_HH

// ## `TypeTag` and compile-time settings (`properties.hh`)
//
// This file defines the `TypeTag` used for the biomineralization example,
// for which we then specialize `properties` to the needs of this setup.
//
// [[content]]
//
// ### Include files
// <details>
//
// This type tag specializes most of the `properties` required for two phase flow with
// multiple components including mineralisation simulations (2pncmin) in DuMuX
// We will use this in the following to inherit the respective properties and
// subsequently specialize those `properties` for our `TypeTag`, which we want to
// modify or for which no meaningful default can be set.
#include <dumux/common/properties.hh>
#include <dumux/material/components/simpleco2.hh>
#include <dumux/porousmediumflow/2pncmin/model.hh>

// We want to use `YaspGrid`, an implementation of the dune grid interface for structured grids:
#include <dune/grid/yaspgrid.hh>

// In this example, we want to discretize the equations with the box scheme
#include <dumux/discretization/box.hh>

// We include the necessary material files
#include <examples/biomineralization/material/fluidsystems/biominsimplechemistry.hh>
#include <examples/biomineralization/material/solidsystems/biominsolids.hh>

// We include the problem and spatial parameters headers used for this simulation.
#include "problem.hh"
#include "spatialparams.hh"

// </details>
//
// ### Definition of the `TypeTag` used for the biomineralization problem
//
// We define a `TypeTag` for our simulation with the name `MICPColumnSimpleChemistry`
// and inherit the `properties` specialized for the type tags `TwoPNCMin` and `BoxModel` respectively.
// This way, most of the `properties` required for this simulations
// using the box scheme are conveniently specialized for our new `TypeTag`.
// However, some properties depend on user choices and no meaningful default value can be set.
// Those properties will be addressed later in this file.
// [[codeblock]]
namespace Dumux::Properties {

// We create new type tag for our simulation which inherits from the 2pncmin model and the box discretization
namespace TTag {
struct MICPColumnSimpleChemistry { using InheritsFrom = std::tuple<TwoPNCMin, BoxModel>; };
} // end namespace TTag
// [[/codeblock]]

// ### Property specializations
//
// In the following piece of code, mandatory `properties` for which no meaningful
// default can be set, are specialized for our type tag `MICPColumnSimpleChemistry`.

// [[codeblock]]
// We set the grid to a 1D Yasp Grid
template<class TypeTag>
struct Grid<TypeTag, TTag::MICPColumnSimpleChemistry>
{ using type = Dune::YaspGrid<1>; };

// We set the problem  used for our simulation, defining boundary and initial conditions (see below)
template<class TypeTag>
struct Problem<TypeTag, TTag::MICPColumnSimpleChemistry>
{ using type = MICPColumnProblemSimpleChemistry<TypeTag>; };

// We set the fluidSystem  used for our simulation
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::MICPColumnSimpleChemistry>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using CO2Impl = Components::SimpleCO2<Scalar>;
    using H2OTabulated = Components::TabulatedComponent<Components::H2O<Scalar>>;
    using type = Dumux::FluidSystems::BioMinSimpleChemistryFluid<Scalar, CO2Impl, H2OTabulated>;
};

// We set the solidSystem  used for our simulation
template<class TypeTag>
struct SolidSystem<TypeTag, TTag::MICPColumnSimpleChemistry>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = SolidSystems::BioMinSolidPhase<Scalar>;
};

// We define the spatial parameters for our simulation. The values are specified in the corresponding spatialparameters header file, which is included above.
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::MICPColumnSimpleChemistry>
{
    using type = ICPSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                  GetPropType<TypeTag, Properties::Scalar>>;
};

// We set the two-phase primary variable formulation used for our simulation
template<class TypeTag>
struct Formulation<TypeTag, TTag::MICPColumnSimpleChemistry>
{ static constexpr auto value = TwoPFormulation::p0s1; };

}// We leave the namespace Dumux::Properties.
// [[/codeblock]]
// [[/content]]
#endif
