// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup RichardsTests
 * \brief The properties of the one-dimensional infiltration problem with a smooth, given solution.
 */
#ifndef DUMUX_RICHARDS_ANALYTICALPROPERTIES_HH
#define DUMUX_RICHARDS_ANALYTICALPROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>

#include <dumux/porousmediumflow/richards/model.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/air.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidsystems/1pgas.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>

#include "spatialparams.hh"
#include "problem.hh"

//////////
// Specify the properties for the analytical problem
//////////
namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct RichardsAnalytical { using InheritsFrom = std::tuple<Richards>; };
struct RichardsAnalyticalCC { using InheritsFrom = std::tuple<RichardsAnalytical, CCTpfaModel>; };
} // end namespace TTag

// Use 2d YaspGrid
template<class TypeTag>
struct Grid<TypeTag, TTag::RichardsAnalytical> { using type = Dune::YaspGrid<2>; };

// Set the physical problem to be solved
template<class TypeTag>
struct Problem<TypeTag, TTag::RichardsAnalytical> { using type = RichardsAnalyticalProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::RichardsAnalytical>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = RichardsAnalyticalSpatialParams<GridGeometry, Scalar>;
};

// Set the fluid system
// -> this is mainly to test that the Richards model also works with
// the immiscible two-phase fluid system. The default fluid system (H2OAir)
// would also work (and the reference solution was created with it)
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::RichardsAnalytical>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using L = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
    using G = FluidSystems::OnePGas<Scalar, Components::Air<Scalar> >;
    using type = FluidSystems::TwoPImmiscible<Scalar, L, G>;
};

// TODO: remove after release (3.6)
// Set the primary variables type
template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::RichardsAnalytical>
{ using type = Dune::FieldVector<GetPropType<TypeTag, Properties::Scalar>, GetPropType<TypeTag, Properties::ModelTraits>::numEq()>; };

} // end namespace Dumux::Properties

#endif
