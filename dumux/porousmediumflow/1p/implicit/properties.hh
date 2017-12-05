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
 * \ingroup OnePModel
 * \file
 *
 * \brief Defines the properties required for the one-phase fully implicit model.
 */
#ifndef DUMUX_1P_PROPERTIES_HH
#define DUMUX_1P_PROPERTIES_HH

#include <dumux/common/properties.hh>

#include <dumux/material/components/nullcomponent.hh>
#include <dumux/material/fluidmatrixinteractions/1p/thermalconductivityaverage.hh>
#include <dumux/material/fluidstates/immiscible.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/fluidsystems/1p.hh>

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/immiscible/localresidual.hh>
#include <dumux/porousmediumflow/nonisothermal/implicit/properties.hh>

#include "indices.hh"
#include "volumevariables.hh"
#include "vtkoutputfields.hh"

namespace Dumux
{
namespace Properties {
//! The type tags for the isothermal & non-isothermal single phase model
NEW_TYPE_TAG(OneP, INHERITS_FROM(PorousMediumFlow));
NEW_TYPE_TAG(OnePNI, INHERITS_FROM(OneP, NonIsothermal));

///////////////////////////////////////////////////////////////////////////
// properties for the isothermal single phase model
///////////////////////////////////////////////////////////////////////////
SET_INT_PROP(OneP, NumEq, 1);                                         //! set the number of equations to 1
SET_INT_PROP(OneP, NumPhases, 1);                                     //! The number of phases in the 1p model is 1
SET_INT_PROP(OneP, NumComponents, 1);                                 //! The number of components in the 1p model is 1
SET_TYPE_PROP(OneP, LocalResidual, ImmiscibleLocalResidual<TypeTag>); //! The local residual function
SET_TYPE_PROP(OneP, VolumeVariables, OnePVolumeVariables<TypeTag>);   //! the VolumeVariables property
SET_BOOL_PROP(OneP, EnableAdvection, true);                           //! The one-phase model considers advection
SET_BOOL_PROP(OneP, EnableMolecularDiffusion, false);                 //! The one-phase model has no molecular diffusion
SET_BOOL_PROP(OneP, EnableEnergyBalance, false);                      //! Isothermal model by default
SET_TYPE_PROP(OneP, Indices, OnePIndices);                            //! The indices required by the isothermal single-phase model
SET_TYPE_PROP(OneP, VtkOutputFields, OnePVtkOutputFields<TypeTag>);   //! Set the vtk output fields specific to this model

//! The single-phase fluid system is used by default
SET_TYPE_PROP(OneP,
              FluidSystem,
              FluidSystems::OneP<typename GET_PROP_TYPE(TypeTag, Scalar), typename GET_PROP_TYPE(TypeTag, Fluid)>);

//! We set a fluid that only throws exceptions.
//! This hopefully makes the user set this property correctly
SET_PROP(OneP, Fluid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, NullComponent<Scalar> > type;
};

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
SET_PROP(OneP, FluidState)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
public:
    typedef ImmiscibleFluidState<Scalar, FluidSystem> type;
};

///////////////////////////////////////////////////////////////////////////
// properties for the non-isothermal single phase model
///////////////////////////////////////////////////////////////////////////
SET_INT_PROP(OnePNI, IsothermalNumEq, 1);                                           //! set number of equations of isothermal model
SET_BOOL_PROP(OnePNI, EnableEnergyBalance, true);                                   //! we do solve for the energy balance here
SET_TYPE_PROP(OnePNI, IsothermalVtkOutputFields, OnePVtkOutputFields<TypeTag>);     //! the isothermal vtk output fields
SET_TYPE_PROP(OnePNI, IsothermalVolumeVariables, OnePVolumeVariables<TypeTag>);     //! Vol vars of the isothermal model
SET_TYPE_PROP(OnePNI, IsothermalLocalResidual, ImmiscibleLocalResidual<TypeTag>);   //! Local residual of the isothermal model
SET_TYPE_PROP(OnePNI, IsothermalIndices, OnePIndices);                              //! Indices of the isothermal model
SET_TYPE_PROP(OnePNI,
              ThermalConductivityModel,
              ThermalConductivityAverage<typename GET_PROP_TYPE(TypeTag, Scalar)>); //! Use the average for effective conductivities

} // end namespace Properties
} // end namespace Dumux

#endif
