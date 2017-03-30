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
#ifndef DUMUX_NAVIER_STOKES_NI_PROPERTY_DEFAULTS_HH
#define DUMUX_NAVIER_STOKES_NI_PROPERTY_DEFAULTS_HH

// TODO clean-up
#include "properties.hh"

#include "model.hh"
#include "volumevariables.hh"
#include "indices.hh"
#include "localresidual.hh"
#include "fluxvariables.hh"
#include "../staggered/problem.hh"
#include "../staggered/propertydefaults.hh"

#include <dumux/implicit/staggered/localresidual.hh>
#include <dumux/material/fluidsystems/gasphase.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/nullcomponent.hh>
#include <dumux/material/fluidsystems/1p.hh>

#include <dumux/material/fluidstates/immiscible.hh>


namespace Dumux
{

namespace Properties
{
// forward declaration
NEW_PROP_TAG(FluxVariables);
NEW_PROP_TAG(FluxVariablesCache);
}
// \{

///////////////////////////////////////////////////////////////////////////
// default property values for the non-isothermal single phase model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

SET_INT_PROP(NavierStokesNI, NumEqCellCenter, 2);

//! the VolumeVariables property
SET_TYPE_PROP(NavierStokesNI, VolumeVariables, NavierStokesNIVolumeVariables<TypeTag>);
SET_TYPE_PROP(NavierStokesNI, Model, NavierStokesNIModel<TypeTag>);
SET_TYPE_PROP(NavierStokesNI, Indices, NavierStokesNIIndices<TypeTag>);

SET_BOOL_PROP(NavierStokesNI, EnableEnergyBalanceStokes, true);

SET_BOOL_PROP(NavierStokesNI, UseMoles, true);

SET_TYPE_PROP(NavierStokesNI, HeatConductionType, FouriersLaw<TypeTag>);

SET_INT_PROP(NavierStokesNI, PhaseIdx, 0); //!< Defines the phaseIdx


////! average is used as default model to compute the effective thermal heat conductivity
//SET_PROP(NavierStokesNI, ThermalConductivityModel)
//{ private :
//    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
//  public:
//    typedef ThermalConductivityAverage<Scalar> type;
//};


/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
//SET_PROP(NavierStokesNI, FluidState)
//{
//    private:
//        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
//        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
//    public:
//        typedef ImmiscibleFluidState<Scalar, FluidSystem> type;
//};

// //! Enable advection
// SET_BOOL_PROP(NavierStokes, EnableAdvection, true);
//
// //! The one-phase model has no molecular diffusion
// SET_BOOL_PROP(NavierStokes, EnableMolecularDiffusion, false);
//
//! Non-Isothermal model by default
//SET_BOOL_PROP(NavierStokesNI, EnableEnergyBalance, true);
//
// //! The indices required by the isothermal single-phase model
// SET_TYPE_PROP(NavierStokes, Indices, NavierStokesCommonIndices<TypeTag>);
//
// //! The weight of the upwind control volume when calculating
// //! fluxes. Use central differences by default.
// SET_SCALAR_PROP(NavierStokes, ImplicitMassUpwindWeight, 0.5);
//
// //! weight for the upwind mobility in the velocity calculation
// //! fluxes. Use central differences by default.
// SET_SCALAR_PROP(NavierStokes, ImplicitMobilityUpwindWeight, 0.5);
//
// //! The fluid system to use by default
// SET_TYPE_PROP(NavierStokes, FluidSystem, Dumux::FluidSystems::OneP<typename GET_PROP_TYPE(TypeTag, Scalar), typename GET_PROP_TYPE(TypeTag, Fluid)>);
//
// SET_PROP(NavierStokes, Fluid)
// { private:
//     typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
// public:
//     typedef FluidSystems::LiquidPhase<Scalar, Dumux::NullComponent<Scalar> > type;
// };
//
// /*!
//  * \brief The fluid state which is used by the volume variables to
//  *        store the thermodynamic state. This should be chosen
//  *        appropriately for the model ((non-)isothermal, equilibrium, ...).
//  *        This can be done in the problem.
//  */
// SET_PROP(NavierStokes, FluidState){
//     private:
//         typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
//         typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
//     public:
//         typedef Dumux::ImmiscibleFluidState<Scalar, FluidSystem> type;
// };
//
// // disable velocity output by default
// SET_BOOL_PROP(NavierStokes, VtkAddVelocity, true);
//
// // enable gravity by default
// SET_BOOL_PROP(NavierStokes, ProblemEnableGravity, true);
//
// SET_BOOL_PROP(NavierStokes, EnableInertiaTerms, true);
//
// SET_BOOL_PROP(NavierStokes, EnableEnergyTransport, true);
//

//////////////////////////////////////////////////////////////////
// Property values for isothermal model required for the general non-isothermal model
//////////////////////////////////////////////////////////////////

// set isothermal Model
// SET_TYPE_PROP(NavierStokesNI, IsothermalModel, NavierStokesModel<TypeTag>);

//set isothermal VolumeVariables
// SET_TYPE_PROP(NavierStokesNI, IsothermalVolumeVariables, NavierStokesVolumeVariables<TypeTag>);

//set isothermal LocalResidual
// SET_TYPE_PROP(NavierStokesNI, IsothermalLocalResidual, ImmiscibleLocalResidual<TypeTag>);

//set isothermal Indices
// SET_TYPE_PROP(NavierStokesNI, IsothermalIndices, NavierStokesCommonIndices<TypeTag>);

//set isothermal NumEq
// SET_INT_PROP(NavierStokesNI, IsothermalNumEq, 1);

//set non-isothermal NumEq
// SET_INT_PROP(NavierStokesNI, NonIsothermalNumEq, 1);


// \}
} // end namespace Properties

} // end namespace Dumux

#endif
