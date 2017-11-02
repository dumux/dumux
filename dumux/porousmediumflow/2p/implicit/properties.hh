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
 * \ingroup TwoPModel
 */
/*!
 * \file
 *
 * \brief Defines the properties required for the two-phase fully implicit model.
 */

#ifndef DUMUX_2P_PROPERTIES_HH
#define DUMUX_2P_PROPERTIES_HH

#include <dumux/common/basicproperties.hh>
#include <dumux/linear/linearsolverproperties.hh>

#include <dumux/material/components/nullcomponent.hh>
#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivitysomerton.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>
#include <dumux/material/fluidstates/immiscible.hh>
#include <dumux/material/spatialparams/implicit.hh>

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/immiscible/localresidual.hh>
#include <dumux/porousmediumflow/nonisothermal/implicit/properties.hh>

#include "indices.hh"
#include "volumevariables.hh"
#include "vtkoutputfields.hh"

namespace Dumux
{

////////////////////////////////
// properties
////////////////////////////////
namespace Properties
{

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tags for the isothermal & non-isothermal two-phase model
NEW_TYPE_TAG(TwoP, INHERITS_FROM(PorousMediumFlow, NumericModel, LinearSolverTypeTag));
NEW_TYPE_TAG(TwoPNI, INHERITS_FROM(TwoP, NonIsothermal));

///////////////////////////////////////////////////////////////////////////
// properties for the isothermal two-phase model
///////////////////////////////////////////////////////////////////////////
SET_INT_PROP(TwoP, NumEq, 2);                                                 //! Set the number of equations to 2
SET_INT_PROP(TwoP, NumPhases, 2);                                             //! The number of phases in the 2p model is 2
SET_INT_PROP(TwoP, NumComponents, 2);                                         //! The number of components in the 2p model is 2
SET_INT_PROP(TwoP, Formulation, TwoPFormulation::pwsn);                       //! Set the default formulation to pWsN
SET_BOOL_PROP(TwoP, EnableAdvection, true);                                   //! Enable advection
SET_BOOL_PROP(TwoP, EnableMolecularDiffusion, false);                         //! The two-phase model has no molecular diffusion
SET_BOOL_PROP(TwoP, EnableEnergyBalance, false);                              //! Isothermal model (non-isothermal type tag is below)
SET_TYPE_PROP(TwoP, LocalResidual, ImmiscibleLocalResidual<TypeTag>);         //! Use the immiscible local residual operator for the 2p model
SET_TYPE_PROP(TwoP, VolumeVariables, TwoPVolumeVariables<TypeTag>);           //! the VolumeVariables property
SET_TYPE_PROP(TwoP, SpatialParams, ImplicitSpatialParams<TypeTag>);           //! The spatial parameters. Use ImplicitSpatialParams by default.
SET_TYPE_PROP(TwoP, VtkOutputFields, TwoPVtkOutputFields<TypeTag>);           //! Set the vtk output fields specific to the twop model

SET_TYPE_PROP(TwoP,
              MaterialLawParams,
              typename GET_PROP_TYPE(TypeTag, MaterialLaw)::Params);          //! Extract material law params from the law itself

SET_TYPE_PROP(TwoP,
              Indices,
              TwoPIndices<TypeTag, GET_PROP_VALUE(TypeTag, Formulation), 0>); //! The indices required by the isothermal 2p model

//! By default, we set a null component as wetting phase
SET_PROP(TwoP, WettingPhase)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, NullComponent<Scalar> > type;
};

//! By default, we set a null component as non-wetting phase
SET_PROP(TwoP, NonwettingPhase)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, NullComponent<Scalar> > type;
};

//! The two-phase model uses the immiscible fluid system
SET_PROP(TwoP, FluidSystem)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, NonwettingPhase) NonwettingPhase;

public:
    typedef FluidSystems::TwoPImmiscible<Scalar,
                                                WettingPhase,
                                                NonwettingPhase> type;
};

//! The two-phase model uses the immiscible fluid state
SET_PROP(TwoP, FluidState)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
public:
    typedef ImmiscibleFluidState<Scalar, FluidSystem> type;
};

////////////////////////////////////////////////////////
// properties for the non-isothermal single phase model
////////////////////////////////////////////////////////
SET_BOOL_PROP(TwoPNI, NiOutputLevel, 0);                                          //! temperature is already written by the isothermal model
SET_INT_PROP(TwoPNI, IsothermalNumEq, 2);                                         //! set isothermal NumEq
SET_TYPE_PROP(TwoPNI, IsothermalVolumeVariables, TwoPVolumeVariables<TypeTag>);   //! set isothermal VolumeVariables
SET_TYPE_PROP(TwoPNI, IsothermalLocalResidual, ImmiscibleLocalResidual<TypeTag>); //! set isothermal LocalResidual

//! set isothermal Indices
SET_PROP(TwoPNI, IsothermalIndices)
{
private:
    enum { Formulation = GET_PROP_VALUE(TypeTag, Formulation) };
public:
    typedef TwoPIndices<TypeTag, Formulation, 0> type;
};

//! Somerton is used as default model to compute the effective thermal heat conductivity
SET_PROP(TwoPNI, ThermalConductivityModel)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
public:
    typedef ThermalConductivitySomerton<Scalar, Indices> type;
};




}
}

#endif
