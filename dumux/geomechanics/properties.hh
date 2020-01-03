// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Geomechanics
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
