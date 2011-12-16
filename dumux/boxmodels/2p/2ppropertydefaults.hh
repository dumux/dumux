/*****************************************************************************
 *   Copyright (C) 2008-2010 by Andreas Lauser                               *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \ingroup TwoPBoxModel
 * \ingroup BoxProperties
 * \ingroup Properties
 * \file
 *
 * \brief Defines default values for the properties required by the
 *        twophase box model.
 */
#ifndef DUMUX_2P_PROPERTY_DEFAULTS_HH
#define DUMUX_2P_PROPERTY_DEFAULTS_HH

#include "2pmodel.hh"
#include "2pproblem.hh"
#include "2pindices.hh"
#include "2pfluxvariables.hh"
#include "2pvolumevariables.hh"
#include "2pproperties.hh"

#include <dumux/material/fluidsystems/gasphase.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/nullcomponent.hh>

#include <dumux/material/MpNcfluidsystems/2pimmisciblefluidsystem.hh>
#include <dumux/material/MpNcfluidstates/immisciblefluidstate.hh>

namespace Dumux
{
namespace Properties
{
//////////////////////////////////////////////////////////////////
// Property defaults
//////////////////////////////////////////////////////////////////
SET_INT_PROP(BoxTwoP, NumEq, 2); //!< set the number of equations to 2
SET_INT_PROP(BoxTwoP, NumPhases, 2); //!< The number of phases in the 2p model is 2

//! Set the default formulation to pWsN
SET_INT_PROP(BoxTwoP,
             Formulation,
             TwoPCommonIndices::pwSn);

//! Use the 2p local jacobian operator for the 2p model
SET_TYPE_PROP(BoxTwoP,
              LocalResidual,
              TwoPLocalResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(BoxTwoP, Model, TwoPModel<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(BoxTwoP, VolumeVariables, TwoPVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(BoxTwoP, FluxVariables, TwoPFluxVariables<TypeTag>);

//! the upwind weight for the mass conservation equations.
SET_SCALAR_PROP(BoxTwoP, MassUpwindWeight, 1.0);

//! The indices required by the isothermal 2p model
SET_TYPE_PROP(BoxTwoP,
              TwoPIndices,
              TwoPIndices<GET_PROP_VALUE(TypeTag, PTAG(Formulation)), 0>);

/*!
 * \brief Set the property for the material parameters by extracting
 *        it from the material law.
 */
SET_TYPE_PROP(BoxTwoP,
              MaterialLawParams,
              typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLaw))::Params);

SET_PROP(BoxTwoP, WettingPhase)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::NullComponent<Scalar> > type;
};

SET_PROP(BoxTwoP, NonWettingPhase)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::NullComponent<Scalar> > type;
};

SET_PROP(BoxTwoP, FluidSystem)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(WettingPhase)) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NonWettingPhase)) NonWettingPhase;

public:
    typedef Dumux::TwoPImmiscibleFluidSystem<Scalar,
                                             WettingPhase,
                                             NonWettingPhase> type;
};

SET_PROP(BoxTwoP, FluidState)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
public:
    typedef ImmiscibleFluidState<Scalar, FluidSystem> type;
};

// enable jacobian matrix recycling by default
SET_BOOL_PROP(BoxTwoP, EnableJacobianRecycling, true);
// enable partial reassembling by default
SET_BOOL_PROP(BoxTwoP, EnablePartialReassemble, true);
// disable velocity output by default
SET_BOOL_PROP(BoxTwoP, EnableVelocityOutput, false);

}
//

}

#endif
