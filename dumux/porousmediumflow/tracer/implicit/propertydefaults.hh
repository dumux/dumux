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
 * \ingroup ImplicitProperties
 * \ingroup TracerModel
 * \file
 *
 * \brief Defines some default values for the properties of the fully implicit tracer model.
 */

#ifndef DUMUX_TRACER_PROPERTY_DEFAULTS_HH
#define DUMUX_TRACER_PROPERTY_DEFAULTS_HH

#include "properties.hh"
#include "model.hh"
#include "volumevariables.hh"
#include "indices.hh"
#include "localresidual.hh"
#include "vtkoutputfields.hh"

#include <dumux/material/spatialparams/implicit1p.hh>
#include <dumux/discretization/stationaryvelocityfield.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivityconstant.hh>

namespace Dumux
{
// \{
namespace Properties
{
//////////////////////////////////////////////////////////////////
// Property values
//////////////////////////////////////////////////////////////////
SET_INT_PROP(Tracer, NumPhases, 1); //!< The number of phases

SET_PROP(Tracer, NumEq) //!< set the number of equations
{
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static const int value = FluidSystem::numComponents;
};

SET_PROP(Tracer, NumComponents) //!< set the number of components
{
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static const int value = FluidSystem::numComponents;
};

SET_BOOL_PROP(Tracer, UseMoles, true); //!< Define that mole fractions are used in the balance equations

//! Use the tracer local residual function for the tracer model
SET_TYPE_PROP(Tracer, LocalResidual, TracerLocalResidual<TypeTag>);

//! define the model
SET_TYPE_PROP(Tracer, Model, TracerModel<TypeTag>);

//! Set the default vtk output fields
SET_TYPE_PROP(Tracer, VtkOutputFields, TracerVtkOutputFields<TypeTag>);

//! define the VolumeVariables
SET_TYPE_PROP(Tracer, VolumeVariables, TracerVolumeVariables<TypeTag>);

//! We use darcy's law as the default for the advective fluxes
SET_TYPE_PROP(Tracer, AdvectionType, StationaryVelocityField<TypeTag>);

//! set gravity flag for compatibility
SET_BOOL_PROP(Tracer, ProblemEnableGravity, false);

//! Set the indices used by the tracer model
SET_TYPE_PROP(Tracer, Indices, TracerIndices<TypeTag>);

//! The spatial parameters to be employed.
//! Use ImplicitSpatialParamsOneP by default.
SET_TYPE_PROP(Tracer, SpatialParams, ImplicitSpatialParamsOneP<TypeTag>);

//! Use simple model with constant tortuosity as pm diffusivity model
SET_PROP(Tracer, EffectiveDiffusivityModel)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = DiffusivityConstant<TypeTag>;
};

// physical processes to be considered by the isothermal model
SET_BOOL_PROP(Tracer, EnableAdvection, true);
SET_BOOL_PROP(Tracer, EnableMolecularDiffusion, true);
SET_BOOL_PROP(Tracer, EnableEnergyBalance, false);

//! default value for the forchheimer coefficient
// Source: Ward, J.C. 1964 Turbulent flow in porous media. ASCE J. Hydraul. Div 90.
//        Actually the Forchheimer coefficient is also a function of the dimensions of the
//        porous medium. Taking it as a constant is only a first approximation
//        (Nield, Bejan, Convection in porous media, 2006, p. 10)
SET_SCALAR_PROP(Tracer, SpatialParamsForchCoeff, 0.55);

/*!
 * \brief default value for tortuosity value (tau) used in macroscopic diffusion
 *
 * Value is 0.5 according to Carman 1937: <i>Fluid flow through granular beds</i>
 * \cite carman1937
 */
SET_SCALAR_PROP(Tracer, TauTortuosity, 0.5);

} // end namespace Properties
// \}
} // end namespace Dumux

#endif
