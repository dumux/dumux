// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PorousmediumflowModels
 * \brief Defines a type tag and some properties for models using the box scheme.
 */

#ifndef DUMUX_POROUSMEDIUM_FLOW_PROPERTIES_HH
#define DUMUX_POROUSMEDIUM_FLOW_PROPERTIES_HH

#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>
#include <dumux/io/vtkoutputmodule.hh>

#include <dumux/porousmediumflow/fluxvariables.hh>
#include <dumux/porousmediumflow/fluxvariablescache.hh>
#include <dumux/porousmediumflow/fluxvariablescachefiller.hh>
#include <dumux/porousmediumflow/nonisothermal/localresidual.hh>
#include <dumux/porousmediumflow/velocityoutput.hh>

#include <dumux/flux/darcyslaw.hh>
#include <dumux/flux/fickslaw.hh>
#include <dumux/flux/fourierslaw.hh>
#include <dumux/flux/dispersionflux.hh>

#include <dumux/material/solidstates/inertsolidstate.hh>
#include <dumux/material/solidsystems/1csolid.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidmatrixinteractions/dispersiontensors/scheidegger.hh>

namespace Dumux {
namespace Properties {

//! Type tag for models involving flow in porous media
// Create new type tags
namespace TTag {
struct PorousMediumFlow { using InheritsFrom = std::tuple<ModelProperties>; };
} // end namespace TTag

//! The flux variables for models involving flow in porous media
template<class TypeTag>
struct FluxVariables<TypeTag, TTag::PorousMediumFlow> { using type = PorousMediumFluxVariables<TypeTag>; };

//! The flux variables cache class for models involving flow in porous media
template<class TypeTag>
struct FluxVariablesCache<TypeTag, TTag::PorousMediumFlow> { using type = PorousMediumFluxVariablesCache<TypeTag>; };

//! The flux variables cache filler (FluxVariablesCache is the data type,
//! the filler knows how to build up the caches for the stencil efficiently)
template<class TypeTag>
struct FluxVariablesCacheFiller<TypeTag, TTag::PorousMediumFlow> { using type = PorousMediumFluxVariablesCacheFiller<TypeTag>; };

//! By default, we use darcy's law for the advective fluxes
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::PorousMediumFlow> { using type = DarcysLaw<TypeTag>; };

//! By default, we use fick's law for the diffusive fluxes
template<class TypeTag>
struct MolecularDiffusionType<TypeTag, TTag::PorousMediumFlow> { using type = FicksLaw<TypeTag>; };

//! Per default, we do not include compositional dispersion effects
template<class TypeTag>
struct EnableCompositionalDispersion<TypeTag, TTag::PorousMediumFlow> { static constexpr bool value = false; };

//! Per default, we do not include thermal dispersion effects
template<class TypeTag>
struct EnableThermalDispersion<TypeTag, TTag::PorousMediumFlow> { static constexpr bool value = false; };

//! By default, we use a diffusive flux for the dispersive fluxes
template<class TypeTag>
struct DispersionFluxType<TypeTag, TTag::PorousMediumFlow> { using type = DiffusiveDispersionFlux<TypeTag>; };

//! By default, we use Scheideggers's law for the dispersive tensor calculation
template<class TypeTag>
struct CompositionalDispersionModel<TypeTag, TTag::PorousMediumFlow> { using type = ScheideggersDispersionTensor<TypeTag>; };

//! By default, we use the same dispersion tensor type as the componsitonal dispersion for the thermal disperion tensor
template<class TypeTag>
struct ThermalDispersionModel<TypeTag, TTag::PorousMediumFlow>
{ using type = GetPropType<TypeTag, Properties::CompositionalDispersionModel>; };

//! By default, we use fourier's law as the default for heat conduction fluxes
template<class TypeTag>
struct HeatConductionType<TypeTag, TTag::PorousMediumFlow> { using type = FouriersLaw<TypeTag>; };

//! By default, parameters are solution-dependent
template<class TypeTag>
struct SolutionDependentAdvection<TypeTag, TTag::PorousMediumFlow> { static constexpr bool value = true; };
template<class TypeTag>
struct SolutionDependentMolecularDiffusion<TypeTag, TTag::PorousMediumFlow> { static constexpr bool value = true; };
template<class TypeTag>
struct SolutionDependentHeatConduction<TypeTag, TTag::PorousMediumFlow> { static constexpr bool value = true; };

//! The default implementation of the energy balance equation for flow problems in porous media.
template<class TypeTag>
struct EnergyLocalResidual<TypeTag, TTag::PorousMediumFlow> { using type = Dumux::EnergyLocalResidual<TypeTag> ; };

//! Velocity output
template<class TypeTag>
struct VelocityOutput<TypeTag, TTag::PorousMediumFlow>
{
    using type = PorousMediumFlowVelocityOutput<GetPropType<TypeTag, Properties::GridVariables>,
                                                GetPropType<TypeTag, Properties::FluxVariables>>;
};

template<class TypeTag>
struct EnableThermalNonEquilibrium<TypeTag, TTag::PorousMediumFlow> { static constexpr bool value = false; };

//! Per default, we disable the box interface solver
template<class TypeTag>
struct EnableBoxInterfaceSolver<TypeTag, TTag::PorousMediumFlow> { static constexpr bool value = false; };

//! per default solid state is inert
template<class TypeTag>
struct SolidState<TypeTag, TTag::PorousMediumFlow>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;
public:
    using type = InertSolidState<Scalar, SolidSystem>;
};

// per default the solid system is inert with one constant component
template<class TypeTag>
struct SolidSystem<TypeTag, TTag::PorousMediumFlow>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using InertComponent = Components::Constant<1, Scalar>;
    using type = SolidSystems::InertSolidPhase<Scalar, InertComponent>;
};

} // namespace Properties
} // namespace Dumux

 #endif
