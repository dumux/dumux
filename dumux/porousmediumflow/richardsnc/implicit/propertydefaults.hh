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
 * \ingroup RichardsNCModel
 * \file
 *
 * \brief Defines some default values for the properties of the
 *        Richards, n-component fully implicit model.
 */

#ifndef DUMUX_RICHARDSNC_PROPERTY_DEFAULTS_HH
#define DUMUX_RICHARDSNC_PROPERTY_DEFAULTS_HH

#include "properties.hh"
#include "model.hh"
#include "volumevariables.hh"
#include "indices.hh"

#include <dumux/porousmediumflow/compositional/localresidual.hh>
#include <dumux/porousmediumflow/richards/implicit/newtoncontroller.hh>
#include <dumux/porousmediumflow/nonisothermal/implicit/propertydefaults.hh>
#include <dumux/material/spatialparams/implicit1p.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh>
#include <dumux/material/fluidmatrixinteractions/1p/thermalconductivityaverage.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/liquidphase2c.hh>
#include <dumux/material/fluidstates/compositional.hh>

namespace Dumux
{
// \{
namespace Properties
{
//////////////////////////////////////////////////////////////////
// Property values
//////////////////////////////////////////////////////////////////
SET_PROP(RichardsNC, NumEq)
{
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static const int value = FluidSystem::numComponents;
};

SET_PROP(RichardsNC, NumPhases)
{
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static const int value = FluidSystem::numPhases;
    static_assert(value == 1, "You can only use one-phasic fluid systems with the Richards n-component model!");
};

SET_PROP(RichardsNC, NumComponents)
{
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static const int value = FluidSystem::numComponents;
};

SET_BOOL_PROP(RichardsNC, UseMoles, true); //!< Define that per default mole fractions are used in the balance equations

//! Use the dedicated local residual
SET_TYPE_PROP(RichardsNC, LocalResidual, CompositionalLocalResidual<TypeTag>);

//! We set the replaceCompIdx to 0, i.e. the first equation is substituted with
//! the total mass balance, i.e. the phase balance
SET_INT_PROP(RichardsNC, ReplaceCompEqIdx, 0);

//! define the model
SET_TYPE_PROP(RichardsNC, Model, RichardsNCModel<TypeTag>);

//! define the VolumeVariables
SET_TYPE_PROP(RichardsNC, VolumeVariables, RichardsNCVolumeVariables<TypeTag>);

//! Smarter newton controller
SET_TYPE_PROP(RichardsNC, NewtonController, RichardsNewtonController<TypeTag>);

/*!
 *\brief The fluid system used by the model.
 *
 * By default this uses the liquid phase fluid system with simple H2O.
 */
SET_PROP(RichardsNC, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = FluidSystems::LiquidPhaseTwoC<TypeTag, SimpleH2O<Scalar>, Constant<TypeTag, Scalar>>;
};

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
SET_PROP(RichardsNC, FluidState)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using type = CompositionalFluidState<Scalar, FluidSystem>;
};

/*!
 * \brief Set type of the parameter objects for the material law
 *
 * By default this is just retrieved from the material law.
 */
SET_TYPE_PROP(RichardsNC, MaterialLawParams, typename GET_PROP_TYPE(TypeTag, MaterialLaw)::Params);

//! Set the indices used
SET_TYPE_PROP(RichardsNC, Indices, RichardsNCIndices<TypeTag>);
//! The spatial parameters to be employed.
//! Use ImplicitSpatialParamsOneP by default.
SET_TYPE_PROP(RichardsNC, SpatialParams, ImplicitSpatialParamsOneP<TypeTag>);

//! The model after Millington (1961) is used for the effective diffusivity
SET_PROP(RichardsNC, EffectiveDiffusivityModel)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = DiffusivityMillingtonQuirk<Scalar>;
};

// enable gravity by default
SET_BOOL_PROP(RichardsNC, ProblemEnableGravity, true);

//physical processes to be considered by the isothermal model
SET_BOOL_PROP(RichardsNC, EnableAdvection, true);
SET_BOOL_PROP(RichardsNC, EnableMolecularDiffusion, true);
SET_BOOL_PROP(RichardsNC, EnableEnergyBalance, false);

//! default value for the forchheimer coefficient
// Source: Ward, J.C. 1964 Turbulent flow in porous media. ASCE J. Hydraul. Div 90.
//        Actually the Forchheimer coefficient is also a function of the dimensions of the
//        porous medium. Taking it as a constant is only a first approximation
//        (Nield, Bejan, Convection in porous media, 2006, p. 10)
SET_SCALAR_PROP(RichardsNC, SpatialParamsForchCoeff, 0.55);

/*!
 * \brief default value for tortuosity value (tau) used in macroscopic diffusion
 *
 * Value is 0.5 according to Carman 1937: <i>Fluid flow through granular beds</i>
 * \cite carman1937
 */
SET_SCALAR_PROP(RichardsNC, TauTortuosity, 0.5);

//! average is used as default model to compute the effective thermal heat conductivity
SET_PROP(RichardsNCNI, ThermalConductivityModel)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = ThermalConductivityAverage<Scalar>;
};

//////////////////////////////////////////////////////////////////
// Property values for isothermal model required for the general non-isothermal model
//////////////////////////////////////////////////////////////////

// set isothermal Model
SET_TYPE_PROP(RichardsNCNI, IsothermalModel, RichardsNCModel<TypeTag>);

//set isothermal VolumeVariables
SET_TYPE_PROP(RichardsNCNI, IsothermalVolumeVariables, RichardsNCVolumeVariables<TypeTag>);

//set isothermal LocalResidual
SET_TYPE_PROP(RichardsNCNI, IsothermalLocalResidual, CompositionalLocalResidual<TypeTag>);

//set isothermal Indices
SET_TYPE_PROP(RichardsNCNI, IsothermalIndices, RichardsNCIndices<TypeTag>);

//set isothermal NumEq
SET_PROP(RichardsNCNI, IsothermalNumEq)
{
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static const int value = FluidSystem::numComponents;
};

} // end namespace Properties
// \}
} // end namespace Dumux

#endif
