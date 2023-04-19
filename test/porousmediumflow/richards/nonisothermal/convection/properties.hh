// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup RichardsTests
 * \brief The properties of the test for the RichardsModel in combination with the NI model for a
 * convection problem.
 */
#ifndef DUMUX_RICHARDS_CONVECTION_PROPERTIES_HH
#define DUMUX_RICHARDS_CONVECTION_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/richards/model.hh>
#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivity/somerton.hh>
#include <dumux/material/fluidsystems/h2on2.hh>

#include "../spatialparams.hh"
#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct RichardsNIConvection { using InheritsFrom = std::tuple<RichardsNI>; };
struct RichardsNIConvectionBox { using InheritsFrom = std::tuple<RichardsNIConvection, BoxModel>; };
struct RichardsNIConvectionCC { using InheritsFrom = std::tuple<RichardsNIConvection, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::RichardsNIConvection> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::RichardsNIConvection> { using type = RichardsNIConvectionProblem<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::RichardsNIConvection> { using type = FluidSystems::H2ON2<GetPropType<TypeTag, Properties::Scalar>, FluidSystems::H2ON2DefaultPolicy</*fastButSimplifiedRelations=*/true>>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::RichardsNIConvection>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = RichardsNISpatialParams<GridGeometry, Scalar>;
};

// TODO: remove after release (3.6)
// Set the primary variables type
template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::RichardsNIConvection>
{ using type = Dune::FieldVector<GetPropType<TypeTag, Properties::Scalar>, GetPropType<TypeTag, Properties::ModelTraits>::numEq()>; };

} // end namespace Dumux::Properties

#endif
