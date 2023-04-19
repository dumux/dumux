// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkModels
 * \brief Defines common properties required for all pore-network models.
 */
#ifndef DUMUX_PNM_COMMON_PROPERTIES_HH
#define DUMUX_PNM_COMMON_PROPERTIES_HH

#include <dumux/common/properties/model.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/porenetwork/gridgeometry.hh>
#include <dumux/flux/porenetwork/fourierslaw.hh>
#include <dumux/porousmediumflow/fluxvariables.hh>

#include <dumux/porenetwork/common/labels.hh>
#include <dumux/porenetwork/common/velocityoutput.hh>

namespace Dumux::Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the pore-network problem
// Create new type tags
namespace TTag {
// TODO Do we want to inherit from box? Or is this a completely new discretization?
struct PoreNetworkModel { using InheritsFrom = std::tuple<ModelProperties, BoxModel>; };
} // end namespace TTag

//////////////////////////////////////////////////////////////////
// New property tags
//////////////////////////////////////////////////////////////////
template<class TypeTag, class MyTypeTag>
struct Labels { using type = UndefinedProperty; }; //!< The pore/throat labels

//////////////////////////////////////////////////////////////////
// Property defaults
//////////////////////////////////////////////////////////////////
//! Set the default for the grid geometry
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::PoreNetworkModel>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = Dumux::PoreNetwork::GridGeometry<Scalar, GridView, enableCache>;
};

template<class TypeTag>
struct HeatConductionType<TypeTag, TTag::PoreNetworkModel> { using type = Dumux::PoreNetwork::PNMFouriersLaw<>; };

//! The labels
template<class TypeTag>
struct Labels<TypeTag, TTag::PoreNetworkModel> { using type = Dumux::PoreNetwork::Labels; };

template<class TypeTag>
struct VelocityOutput<TypeTag, TTag::PoreNetworkModel>
{
private:
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
public:
    using type = Dumux::PoreNetwork::VelocityOutput<GridVariables, FluxVariables>;
};

template<class TypeTag>
struct EnableThermalNonEquilibrium<TypeTag, TTag::PoreNetworkModel> { static constexpr bool value = false; };

} // end Dumux::Properties

#endif
