// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup RichardsTests
 * \brief The properties of the test for the RichardsModel in combination with the NI model for
 * evaporation.
 */
#ifndef DUMUX_RICHARDS_EVAPORATION_PROPERTIES_HH
#define DUMUX_RICHARDS_EVAPORATION_PROPERTIES_HH
/**
 */
#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/richardsextended/model.hh>
#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivity/somerton.hh>
#include <dumux/material/fluidsystems/h2on2.hh>

#include "../spatialparams.hh"
#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct RichardsNIEvaporation { using InheritsFrom = std::tuple<ExtendedRichardsNI>; };
struct RichardsNIEvaporationBox { using InheritsFrom = std::tuple<RichardsNIEvaporation, BoxModel>; };
struct RichardsNIEvaporationCC { using InheritsFrom = std::tuple<RichardsNIEvaporation, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::RichardsNIEvaporation> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::RichardsNIEvaporation> { using type = RichardsNIEvaporationProblem<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::RichardsNIEvaporation> { using type = FluidSystems::H2ON2<GetPropType<TypeTag, Properties::Scalar>, FluidSystems::H2ON2DefaultPolicy</*fastButSimplifiedRelations=*/true>>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::RichardsNIEvaporation>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = RichardsNISpatialParams<GridGeometry, Scalar>;
};

} // end namespace Dumux::Properties

#endif
