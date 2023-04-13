// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPTwoCTests
 * \brief The properties of a evaporation problem where two components with constant properties mix and evaporate.
 */

#ifndef DUMUX_EVAPORATION_CONSTANT_COMPONENT_PROPERTIES_HH
#define DUMUX_EVAPORATION_CONSTANT_COMPONENT_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/2p2c/model.hh>
#include <dumux/material/components/simpleh2o.hh>

#include "spatialparams.hh"
#include "constant2p2cfluidsystem.hh"
#include "problem.hh"

#ifndef WETTINGCOMPONENT
#define WETTINGCOMPONENT Components::Constant<2,Scalar>
#endif

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct EvaporationConstantComponent { using InheritsFrom = std::tuple<TwoPTwoCNI>; };
struct EvaporationConstantComponentBox { using InheritsFrom = std::tuple<EvaporationConstantComponent, BoxModel>; };
struct EvaporationConstantComponentCCTpfa { using InheritsFrom = std::tuple<EvaporationConstantComponent, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::EvaporationConstantComponent> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::EvaporationConstantComponent> { using type = EvaporationConstantComponentProblem<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::EvaporationConstantComponent>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::ConstantTwoPTwoCFluidsystem<Scalar, WETTINGCOMPONENT, Components::Constant<3, Scalar> >;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::EvaporationConstantComponent>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = EvaporationConstantComponentSpatialParams<GridGeometry, Scalar>;
};

// Define whether mole(true) or mass (false) fractions are used
template<class TypeTag>
struct UseMoles<TypeTag, TTag::EvaporationConstantComponent> { static constexpr bool value = true; };

//! Set the  formulation to pw-Sn
template<class TypeTag>
struct Formulation<TypeTag, TTag::EvaporationConstantComponent>
{ static constexpr auto value = TwoPFormulation::p0s1; };

// Enable caching
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::EvaporationConstantComponent> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::EvaporationConstantComponent> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::EvaporationConstantComponent> { static constexpr bool value = true; };

} // end namespace Dumux::Properties

#endif
