// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Default properties for the vector-potential Boussinesq model.
 *
 * Two TypeTags are provided:
 *   BoussinesqVorticityModel     – physics (ModelTraits, VolumeVariables,
 *                                  LocalResidual, IOFields, VelocityOutput)
 *   BoussinesqVorticityModelBox  – adds Box (CVFE) discretisation
 *
 * The model defaults to a single transported component and 2D geometry.
 * To use 3D or multiple components, specialize ModelTraits in your test TypeTag:
 *
 *   template<class TypeTag>
 *   struct ModelTraits<TypeTag, TTag::MyTest>
 *   { using type = BoussinesqVorticityModelTraits<nComp, 3>; };
 */
#ifndef DUMUX_BOUSSINESQ_VORTICITY_MODEL_HH
#define DUMUX_BOUSSINESQ_VORTICITY_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>
#include <dumux/discretization/box.hh>

#include "indices.hh"
#include "modeltraits.hh"
#include "volumevariables.hh"
#include "localresidual.hh"
#include "iofields.hh"
#include "velocityoutput.hh"

namespace Dumux::Properties {

// -------------------------------------------------------------------------
// TypeTags
// -------------------------------------------------------------------------
namespace TTag {

struct BoussinesqVorticityModel
{ using InheritsFrom = std::tuple<ModelProperties>; };

struct BoussinesqVorticityModelBox
{ using InheritsFrom = std::tuple<BoussinesqVorticityModel, BoxModel>; };

// Backward-compat aliases
using BoussinesqStreamfunctionModel    = BoussinesqVorticityModel;
using BoussinesqStreamfunctionModelBox = BoussinesqVorticityModelBox;

} // end namespace TTag

// -------------------------------------------------------------------------
// Model traits — defaults to 1 component, 2D.
// Override in the test-specific TypeTag to change nComp or dim.
// -------------------------------------------------------------------------
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::BoussinesqVorticityModel>
{ using type = BoussinesqVorticityModelTraits<1, 2>; };

// -------------------------------------------------------------------------
// Volume variables
// -------------------------------------------------------------------------
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::BoussinesqVorticityModel>
{
    struct Traits
    {
        using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
        using ModelTraits      = GetPropType<TypeTag, Properties::ModelTraits>;
    };
    using type = BoussinesqVorticityVolumeVariables<Traits>;
};

// -------------------------------------------------------------------------
// Local residual
// -------------------------------------------------------------------------
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::BoussinesqVorticityModel>
{ using type = BoussinesqVorticityLocalResidual<TypeTag>; };

// -------------------------------------------------------------------------
// IO fields
// -------------------------------------------------------------------------
template<class TypeTag>
struct IOFields<TypeTag, TTag::BoussinesqVorticityModel>
{ using type = BoussinesqVorticityIOFields; };

// -------------------------------------------------------------------------
// Velocity output (Box-specific)
// -------------------------------------------------------------------------
template<class TypeTag>
struct VelocityOutput<TypeTag, TTag::BoussinesqVorticityModelBox>
{ using type = BoussinesqVorticityVelocityOutput<GetPropType<TypeTag, Properties::GridVariables>>; };

} // end namespace Dumux::Properties

#endif
