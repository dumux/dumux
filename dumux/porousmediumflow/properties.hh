// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup Properties
 * \file
 *
 * \brief Defines a type tag and some properties for models using the box scheme.
 */

#ifndef DUMUX_POROUSMEDIUM_FLOW_PROPERTIES_HH
#define DUMUX_POROUSMEDIUM_FLOW_PROPERTIES_HH

#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>
#include <dumux/io/vtkoutputmodule.hh>

#include <dumux/porousmediumflow/fluxvariables.hh>
#include <dumux/porousmediumflow/fluxvariablescache.hh>
#include <dumux/porousmediumflow/nonisothermal/localresidual.hh>
#include <dumux/porousmediumflow/compositional/primaryvariableswitch.hh>
#include <dumux/porousmediumflow/velocityoutput.hh>

#include <dumux/discretization/darcyslaw.hh>
#include <dumux/discretization/fickslaw.hh>
#include <dumux/discretization/fourierslaw.hh>

namespace Dumux {
namespace Properties {

//! Type tag for models involving flow in porous media
NEW_TYPE_TAG(PorousMediumFlow, INHERITS_FROM(ModelProperties));

//! The flux variables for models involving flow in porous media
SET_TYPE_PROP(PorousMediumFlow, FluxVariables, PorousMediumFluxVariables<TypeTag>);

//! The flux variables cache class for models involving flow in porous media
SET_TYPE_PROP(PorousMediumFlow, FluxVariablesCache, PorousMediumFluxVariablesCache<TypeTag>);

//! By default, we use darcy's law for the advective fluxes
SET_TYPE_PROP(PorousMediumFlow, AdvectionType, DarcysLaw<TypeTag>);

//! By default, we use fick's law for the diffusive fluxes
SET_TYPE_PROP(PorousMediumFlow, MolecularDiffusionType, FicksLaw<TypeTag>);

//! By default, we use fourier's law as the default for heat conduction fluxes
SET_TYPE_PROP(PorousMediumFlow, HeatConductionType, FouriersLaw<TypeTag>);

//! By default, parameters are solution-dependent
SET_BOOL_PROP(PorousMediumFlow, SolutionDependentAdvection, true);
SET_BOOL_PROP(PorousMediumFlow, SolutionDependentMolecularDiffusion, true);
SET_BOOL_PROP(PorousMediumFlow, SolutionDependentHeatConduction, true);

//! The default implementation of the energy balance equation for flow problems in porous media.
SET_TYPE_PROP(PorousMediumFlow, EnergyLocalResidual, EnergyLocalResidual<TypeTag> );

//! Velocity output
SET_TYPE_PROP(PorousMediumFlow, VelocityOutput, PorousMediumFlowVelocityOutput<TypeTag>);

//! By default, we set an empty primary variables switch
SET_TYPE_PROP(PorousMediumFlow, PrimaryVariableSwitch, NoPrimaryVariableSwitch);

SET_BOOL_PROP(PorousMediumFlow, EnableThermalNonEquilibrium, false);

//! Per default, we disable the box interface solver
SET_BOOL_PROP(PorousMediumFlow, EnableBoxInterfaceSolver, false);
} // namespace Properties
} // namespace Dumux

 #endif
