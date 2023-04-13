// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \brief Defines a type tag and some properties for free flow models.
 */

#ifndef DUMUX_FREE_FLOW_PROPERTIES_HH
#define DUMUX_FREE_FLOW_PROPERTIES_HH

#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>
#include <dumux/flux/fourierslaw.hh>

#include "turbulencemodel.hh"
#include "spatialparams.hh"

namespace Dumux {
namespace Properties {

//! Type tag for free-flow models
// Create new type tags
namespace TTag {
struct FreeFlow { using InheritsFrom = std::tuple<ModelProperties>; };
} // end namespace TTag

//! Use Fourier's Law as default heat conduction type
template<class TypeTag>
struct HeatConductionType<TypeTag, TTag::FreeFlow> { using type = FouriersLaw<TypeTag>; };

// Set the default spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::FreeFlow>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FreeFlowDefaultSpatialParams<GridGeometry, Scalar>;
};

} // namespace Properties
} // namespace Dumux

 #endif
