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
 * \ingroup OnePNCModel
 * \file
 *
 * \brief Defines some default values for the properties required by the
 *        one-phase n-component fully implicit model.
 */
#ifndef DUMUX_1PNC_PROPERTY_DEFAULTS_HH
#define DUMUX_1PNC_PROPERTY_DEFAULTS_HH

#include "indices.hh"
#include "model.hh"
//#include "fluxvariables.hh"
#include "volumevariables.hh"
#include "properties.hh"
#include "newtoncontroller.hh"

#include <dumux/porousmediumflow/compositional/localresidual.hh>
#include <dumux/porousmediumflow/nonisothermal/implicit/propertydefaults.hh>
#include <dumux/material/spatialparams/implicit1p.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh>
#include <dumux/porousmediumflow/implicit/darcyfluxvariables.hh>
#include <dumux/material/fluidmatrixinteractions/1p/thermalconductivityaverage.hh>
#include <dumux/material/fluidstates/compositional.hh>

namespace Dumux
{

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Property values
//////////////////////////////////////////////////////////////////

/*!
 * \brief Set the property for the number of equations: For each existing component one equation has to be solved.
 */
SET_PROP(OnePNC, NumEq)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
public:
    static const int value = FluidSystem::numComponents;
};


SET_INT_PROP(OnePNC, NumPhases,1); //!< The number of phases in the 1pnc model is 1

SET_BOOL_PROP(OnePNC, UseMoles, true); //!< Define that mole fractions are used in the balance equations

/*!
 * \brief Set the property for the number of components.
 *
 * We just forward the number from the fluid system
 *
 */
SET_PROP(OnePNC, NumComponents)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
    static const int value = FluidSystem::numComponents;
};

//! Set as default that no component mass balance is replaced by the total mass balance
SET_PROP(OnePNC, ReplaceCompEqIdx)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
    static const int value = FluidSystem::numComponents;
};


/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
SET_PROP(OnePNC, FluidState){
    private:
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    public:
        typedef Dumux::CompositionalFluidState<Scalar, FluidSystem> type;
};

//! Use the 1pnc local residual
SET_TYPE_PROP(OnePNC, LocalResidual, CompositionalLocalResidual<TypeTag>);

//! Use the 1pnc newton controller
SET_TYPE_PROP(OnePNC, NewtonController, OnePNCNewtonController<TypeTag>);

//! the Model property
SET_TYPE_PROP(OnePNC, Model, OnePNCModel<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(OnePNC, VolumeVariables, OnePNCVolumeVariables<TypeTag>);

//! the FluxVariables property
// SET_TYPE_PROP(OnePNC, FluxVariables, OnePNCFluxVariables<TypeTag>);

//! define the base flux variables to realize Darcy flow
// SET_TYPE_PROP(OnePNC, BaseFluxVariables, ImplicitDarcyFluxVariables<TypeTag>);

//! the upwind weight for the mass conservation equations.
// SET_SCALAR_PROP(OnePNC, ImplicitMassUpwindWeight, 1.0);

// //! Set default mobility upwind weight to 1.0, i.e. fully upwind
// SET_SCALAR_PROP(OnePNC, ImplicitMobilityUpwindWeight, 1.0);

//! The indices required by the isothermal 2pnc model
SET_TYPE_PROP(OnePNC, Indices, OnePNCIndices <TypeTag, /*PVOffset=*/0>);

//! Set the phaseIndex per default to zero (important for one-phase fluidsystems).
SET_INT_PROP(OnePNC, PhaseIdx, 0);

//! Use the ImplicitSpatialParams by default
SET_TYPE_PROP(OnePNC, SpatialParams, ImplicitSpatialParamsOneP<TypeTag>);

//! Enable gravity by default
SET_BOOL_PROP(OnePNC, ProblemEnableGravity, false);

//! Disable velocity output by default
SET_BOOL_PROP(OnePNC, VtkAddVelocity, false);

//! Somerton is used as default model to compute the effective thermal heat conductivity
// SET_PROP(OnePNCNI, ThermalConductivityModel)
// {
// private:
//     typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
//     typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
// public:
//     typedef ThermalConductivityAverage<Scalar> type;
// };

//! average is used as default model to compute the effective thermal heat conductivity
SET_PROP(OnePNC, ThermalConductivityModel)
{ private :
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
  public:
    typedef ThermalConductivityAverage<Scalar> type;
};

//! The model after Millington (1961) is used for the effective diffusivity
SET_PROP(OnePNC, EffectiveDiffusivityModel)
{ private :
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
 public:
    typedef DiffusivityMillingtonQuirk<Scalar> type;
};


//! Enable advection
SET_BOOL_PROP(OnePNC, EnableAdvection, true);

//! Enable molecular diffusion
SET_BOOL_PROP(OnePNC, EnableMolecularDiffusion, true);

//! Isothermal model by default
SET_BOOL_PROP(OnePNC, EnableEnergyBalance, false);
//////////////////////////////////////////////////////////////////
// Property values for isothermal model required for the general non-isothermal model
//////////////////////////////////////////////////////////////////

// set isothermal Model
SET_TYPE_PROP(OnePNCNI, IsothermalModel, OnePNCModel<TypeTag>);

//set isothermal VolumeVariables
SET_TYPE_PROP(OnePNC, IsothermalVolumeVariables, OnePNCVolumeVariables<TypeTag>);

//set isothermal VolumeVariables
SET_TYPE_PROP(OnePNCNI, IsothermalVolumeVariables, OnePNCVolumeVariables<TypeTag>);

//set isothermal LocalResidual
SET_TYPE_PROP(OnePNCNI, IsothermalLocalResidual, CompositionalLocalResidual<TypeTag>);

//set isothermal Indices
SET_TYPE_PROP(OnePNCNI, IsothermalIndices, OnePNCIndices<TypeTag, /*PVOffset=*/0>);

//set isothermal NumEq
SET_PROP(OnePNCNI, IsothermalNumEq)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
    static const int value = FluidSystem::numComponents;
};

//physical processes to be considered by the non-isothermal model
SET_BOOL_PROP(OnePNCNI, EnableAdvection, true);
SET_BOOL_PROP(OnePNCNI, EnableMolecularDiffusion, true);
SET_BOOL_PROP(OnePNCNI, EnableEnergyBalance, true);

}
}

#endif
