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

#include <dumux/material/solidstates/inertsolidstate.hh>
#include <dumux/material/solidsystems/1csolid.hh>
#include <dumux/material/components/constant.hh>

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
