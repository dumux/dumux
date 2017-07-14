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
 * \ingroup TwoPMimeticModel
 * \file
 *
 * \brief Defines the properties required for the one-phase fully implicit model.
 */
#ifndef DUMUX_2P_MIMETIC_PROPERTY_DEFAULTS_HH
#define DUMUX_2P_MIMETIC_PROPERTY_DEFAULTS_HH

#include "properties.hh"

#include "model.hh"
#include "../implicit/volumevariables.hh"
#include "indices.hh"
#include <dumux/porousmediumflow/immiscible/mimetic/localresidual.hh>
#include "problem.hh"
#include <dumux/implicit/staggered/localresidual.hh>
#include <dumux/material/fluidsystems/gasphase.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/nullcomponent.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>
#include <dumux/material/fluidstates/immiscible.hh>
#include <dumux/material/spatialparams/implicit.hh>


namespace Dumux
{

// \{

///////////////////////////////////////////////////////////////////////////
// default property values for the isothermal single phase model
///////////////////////////////////////////////////////////////////////////
namespace Properties {
SET_INT_PROP(TwoPMimetic, NumEqCellCenter, 2); //!< set the number of equations to 2
SET_INT_PROP(TwoPMimetic, NumEqFace, 2); //!< set the number of equations to 2
SET_INT_PROP(TwoPMimetic, NumPhases, 2); //!< The number of phases in the 2p model is 2

//! Set the default formulation to pWsN
SET_INT_PROP(TwoPMimetic,
             Formulation,
             TwoPMimeticFormulation::pwsn);

//! The local residual function
SET_TYPE_PROP(TwoPMimetic, LocalResidual, ImmiscibleMimeticLocalResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(TwoPMimetic, Model, TwoPMimeticModel<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(TwoPMimetic, VolumeVariables, TwoPVolumeVariables<TypeTag>);

//! Enable advection
SET_BOOL_PROP(TwoPMimetic, EnableAdvection, true);

//! The one-phase model has no molecular diffusion
SET_BOOL_PROP(TwoPMimetic, EnableMolecularDiffusion, false);

//! Isothermal model by default
SET_BOOL_PROP(TwoPMimetic, EnableEnergyBalance, false);

//! The indices required by the isothermal 2p model
SET_TYPE_PROP(TwoPMimetic,
              Indices,
              TwoPMimeticIndices<TypeTag, GET_PROP_VALUE(TypeTag, Formulation), 0>);

/*!
 * \brief Set the property for the material parameters by extracting
 *        it from the material law.
 */
SET_TYPE_PROP(TwoPMimetic,
              MaterialLawParams,
              typename GET_PROP_TYPE(TypeTag, MaterialLaw)::Params);

SET_PROP(TwoPMimetic, WettingPhase)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, NullComponent<Scalar> > type;
};

SET_PROP(TwoPMimetic, NonwettingPhase)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, NullComponent<Scalar> > type;
};

SET_PROP(TwoPMimetic, FluidSystem)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, NonwettingPhase) NonwettingPhase;

public:
    typedef FluidSystems::TwoPImmiscible<Scalar,
                                                WettingPhase,
                                                NonwettingPhase> type;
};

SET_PROP(TwoPMimetic, FluidState)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
public:
    typedef ImmiscibleFluidState<Scalar, FluidSystem> type;
};
// disable velocity output by default
SET_BOOL_PROP(TwoPMimetic, VtkAddVelocity, false);

// enable gravity by default
SET_BOOL_PROP(TwoPMimetic, ProblemEnableGravity, true);
;

SET_PROP(TwoPMimetic, BoundaryValues)
{
private:
    using CellCenterBoundaryValues = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FaceBoundaryValues = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);
public:
    using type = StaggeredPrimaryVariables<TypeTag, CellCenterBoundaryValues, FaceBoundaryValues>;
};

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

// \}
} // end namespace Properties

} // end namespace Dumux

#endif
