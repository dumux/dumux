// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup GeomechanicsModels
 * \brief Defines a type tag and some properties for geomechanical DuMuX models.
 */

#ifndef DUMUX_GEOMECHANICS_PROPERTIES_HH
#define DUMUX_GEOMECHANICS_PROPERTIES_HH

#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/solidstates/inertsolidstate.hh>
#include <dumux/material/solidsystems/1csolid.hh>
#include <dumux/flux/hookeslaw.hh>

#include "stressvariablescache.hh"
#include "velocityoutput.hh"

namespace Dumux {
namespace Properties {

//! Type tag for geomechanical models
// Create new type tags
namespace TTag {
struct Geomechanics { using InheritsFrom = std::tuple<ModelProperties>; };
} // end namespace TTag

//! The flux variables cache class for models involving flow in porous media
template<class TypeTag>
struct FluxVariablesCache<TypeTag, TTag::Geomechanics>
{
    using type = StressVariablesCache< GetPropType<TypeTag, Properties::Scalar>,
                                       GetPropType<TypeTag, Properties::GridGeometry> >;
};

//! The (currently empty) velocity output
template<class TypeTag>
struct VelocityOutput<TypeTag, TTag::Geomechanics> { using type = GeomechanicsVelocityOutput<GetPropType<TypeTag, Properties::GridVariables>>; };

//! The solid state must be inert
template<class TypeTag>
struct SolidState<TypeTag, TTag::Geomechanics>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;
public:
    using type = InertSolidState<Scalar, SolidSystem>;
};

//! Per default we use one constant component in the inert solid system
template<class TypeTag>
struct SolidSystem<TypeTag, TTag::Geomechanics>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using InertComponent = Components::Constant<1, Scalar>;
public:
    using type = SolidSystems::InertSolidPhase<Scalar, InertComponent>;
};
} // namespace Properties
} // namespace Dumux

 #endif
