// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPNCMinTests
 * \brief The properties of the problem where brine is evaporating at the top boundary. The system is closed at the remaining boundaries.
 */
#ifndef DUMUX_SALINIZATION_PROPERTIES_HH
#define DUMUX_SALINIZATION_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/porousmediumflow/2pncmin/model.hh>
#include <dumux/material/fluidsystems/brineair.hh>

#include <dumux/material/components/nacl.hh>
#include <dumux/material/components/granite.hh>
#include <dumux/material/solidsystems/compositionalsolidphase.hh>

#include "spatialparams.hh"
#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct Salinization { using InheritsFrom = std::tuple<TwoPNCMinNI>; };
struct SalinizationBox { using InheritsFrom = std::tuple<Salinization, BoxModel>; };
struct SalinizationCCTpfa { using InheritsFrom = std::tuple<Salinization, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::Salinization>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<Scalar, 2>>;
};

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Salinization> { using type = SalinizationProblem<TypeTag>; };

// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Salinization>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::BrineAir<Scalar, Components::H2O<Scalar>>;
};

template<class TypeTag>
struct SolidSystem<TypeTag, TTag::Salinization>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ComponentOne = Components::NaCl<Scalar>;
    using ComponentTwo = Components::Granite<Scalar>;
    static constexpr int numInertComponents = 1;
    using type = SolidSystems::CompositionalSolidPhase<Scalar, ComponentOne, ComponentTwo, numInertComponents>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Salinization>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = SalinizationSpatialParams<GridGeometry, Scalar>;
};

// Set properties here to override the default property settings
template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::Salinization> { static constexpr int value = 1; }; //!< Replace gas balance by total mass balance
template<class TypeTag>
struct Formulation<TypeTag, TTag::Salinization>
{ static constexpr auto value = TwoPFormulation::p0s1; };

// Set constant extrusion factor via property (for Box)
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::SalinizationBox>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    
    struct GGTraits : public BoxDefaultGridGeometryTraits<GridView>
    { using Extrusion = ConstantExtrusion<0.054977871437821>; };
public:
    using type = BoxFVGridGeometry<Scalar, GridView, enableCache, GGTraits>;
};

// Set constant extrusion factor via property (for CCTpfa)
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::SalinizationCCTpfa>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    
    struct GGTraits : public CCTpfaDefaultGridGeometryTraits<GridView>
    { using Extrusion = ConstantExtrusion<0.054977871437821>; };
public:
    using type = CCTpfaFVGridGeometry<GridView, enableCache, GGTraits>;
};

} // end namespace Dumux::Properties

#endif
