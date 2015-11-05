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
 * \ingroup ElOnePTwoCBoxModel
 * \file
 *
 * \brief Defines the properties required for the one-phase two-component
 * linear-elastic model.
 *
 * This class inherits from the properties of the one-phase two-component model and
 * from the properties of the simple linear-elastic model
 */

#ifndef DUMUX_ELASTIC1P2C_PROPERTY_DEFAULTS_HH
#define DUMUX_ELASTIC1P2C_PROPERTY_DEFAULTS_HH

#include "el1p2cproperties.hh"

#include "el1p2cmodel.hh"
#include "el1p2clocalresidual.hh"
#include "el1p2clocaljacobian.hh"
#include "el1p2cfluxvariables.hh"
#include "el1p2celementvolumevariables.hh"
#include "el1p2cvolumevariables.hh"
#include "el1p2cindices.hh"
#include <dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh>
#include <dumux/material/fluidstates/compositionalfluidstate.hh>


namespace Dumux
{
// \{
namespace Properties
{
//////////////////////////////////////////////////////////////////
// Property defaults
//////////////////////////////////////////////////////////////////
//!< set the number of equations to the space dimension of the problem
SET_PROP(BoxElasticOnePTwoC, NumEq)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum{dim = GridView::dimension};
public:
    static const int value = dim + 2;
};

SET_INT_PROP(BoxElasticOnePTwoC, NumPhases, 1); //!< The number of phases in the 1p2c model is 1
SET_INT_PROP(BoxElasticOnePTwoC, NumComponents, 2); //!< The number of components in the 1p2c model is 2

//! Use the linear elasticity local residual function for the elasticity model
SET_TYPE_PROP(BoxElasticOnePTwoC,
              LocalResidual,
              ElOnePTwoCLocalResidual<TypeTag>);

//! Use the linear elasticity local residual function for the elasticity model
SET_TYPE_PROP(BoxElasticOnePTwoC,
              LocalJacobian,
              ElOnePTwoCLocalJacobian<TypeTag>);

//! define the model
SET_TYPE_PROP(BoxElasticOnePTwoC, Model, ElOnePTwoCModel<TypeTag>);

//! define the ElementVolumeVariables
SET_TYPE_PROP(BoxElasticOnePTwoC, ElementVolumeVariables, Dumux::ElOnePTwoCElementVolumeVariables<TypeTag>);

//! define the VolumeVariables
SET_TYPE_PROP(BoxElasticOnePTwoC, VolumeVariables, Dumux::ElOnePTwoCVolumeVariables<TypeTag>);

//! define the FluxVariables
SET_TYPE_PROP(BoxElasticOnePTwoC, FluxVariables, ElOnePTwoCFluxVariables<TypeTag>);

//! Set the indices used by the linear elasticity model
SET_TYPE_PROP(BoxElasticOnePTwoC, Indices, ElOnePTwoCIndices<TypeTag>);

//! Set the phaseIndex per default to zero (important for two-phase fluidsystems).
SET_INT_PROP(BoxElasticOnePTwoC, PhaseIdx, 0);

SET_PROP(BoxElasticOnePTwoC, FluidState){
    private:
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    public:
        typedef Dumux::CompositionalFluidState<Scalar, FluidSystem> type;
};

//! set default upwind weights to 1.0, i.e. fully upwind
SET_SCALAR_PROP(BoxElasticOnePTwoC, ImplicitMassUpwindWeight, 1.0);
SET_SCALAR_PROP(BoxElasticOnePTwoC, ImplicitMobilityUpwindWeight, 1.0);

// enable gravity by default
SET_BOOL_PROP(BoxElasticOnePTwoC, ProblemEnableGravity, true);

// enable gravity by default
SET_BOOL_PROP(BoxElasticOnePTwoC, ImplicitWithStabilization, true);

//! The model after Millington (1961) is used for the effective diffusivity
SET_PROP(BoxElasticOnePTwoC, EffectiveDiffusivityModel)
{ private :
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
 public:
    typedef DiffusivityMillingtonQuirk<Scalar> type;
};

// write the stress and displacement output according to rock mechanics sign convention (compressive stresses > 0)
SET_BOOL_PROP(BoxElasticOnePTwoC, VtkRockMechanicsSignConvention, false);
}
}

#endif

