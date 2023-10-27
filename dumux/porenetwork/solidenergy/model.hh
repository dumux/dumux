// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PNMSolidEnergyModel
 * \brief The energy balance equation for a porous solid in pore-networks based on PorousMediumFlow-SolidEnergy
 */

#ifndef DUMUX_PNM_SOLID_ENERGY_MODEL_HH
#define DUMUX_PNM_SOLID_ENERGY_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/porenetwork/properties.hh>

#include <dumux/porousmediumflow/solidenergy/model.hh>
#include <dumux/flux/porenetwork/grainfourierslaw.hh>

#include "spatialparams.hh"
#include "volumevariables.hh"
#include "fluxvariablescache.hh"
#include "iofields.hh"

/*!
 * \ingroup PNMSolidEnergyModel
 * \brief The energy balance equation for a porous solid
 *
 * The energy balance is described by the following equation:
 \f[
   \frac{ \partial n c_p \varrho T}{\partial t}
   - \nabla \cdot \left\lbrace \lambda_\text{pm} \nabla T \right\rbrace = q,
 \f]
 * where \f$n\f$ is the volume fraction of the conducting material, \f$c_p\f$ its specific heat capacity,
 * \f$\varrho\f$ its density, \f$T\f$ the temperature, and \f$\lambda\f$ the heat conductivity of the porous solid.
*/

///////////////////////////////////////////////////////////////////////////
// properties for the solid-energy model
///////////////////////////////////////////////////////////////////////////
namespace Dumux::Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tags for the implicit solid-energy problems
// Create new type tags
namespace TTag {
struct PNMSolidEnergy{ using InheritsFrom = std::tuple<PoreNetworkModel, SolidEnergy>; };
}

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::PNMSolidEnergy>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;

    using Traits = SolidEnergyVolumeVariablesTraits<PV, SSY, SST, MT>;
public:
    using type = PoreNetwork::SolidEnergyVolumeVariables<Traits>;
};

template<class TypeTag>
struct HeatConductionType<TypeTag, TTag::PNMSolidEnergy>
{ using type = PoreNetwork::TruncatedPyramidGrainFouriersLaw<GetPropType<TypeTag, Properties::Scalar>>; };

//! The spatial parameters to be employed.
//! Use PNMOnePSpatialParams by default.
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::PNMSolidEnergy>
{
private:
    using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = PoreNetwork::SolidEnergySpatialParams<FVGridGeometry, Scalar>;
};

//! The flux variables cache
template<class TypeTag>
struct FluxVariablesCache<TypeTag, TTag::PNMSolidEnergy>
{ using type = PoreNetwork::SolidEnergyFluxVariablesCache<GetPropType<TypeTag, Properties::Scalar>>; };

//! Set the vtk output fields specific to this model
template<class TypeTag>
struct IOFields<TypeTag, TTag::PNMSolidEnergy>
{ using type = PoreNetwork::SolidEnergyIOFields; };

} // namespace Dumux::Properties

#endif
