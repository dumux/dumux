// $Id$
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
 */
/*!
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
#include "2pfluidstate.hh"
#include "2pproperties.hh"
#include <dumux/material/fluidsystems/2p_system.hh>

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

//! the upwind factor for the mobility.
SET_SCALAR_PROP(BoxTwoP, MobilityUpwindAlpha, 1.0);

//! The indices required by the isothermal 2p model
SET_PROP(BoxTwoP, TwoPIndices)
{
    typedef TwoPIndices<GET_PROP_VALUE(TypeTag, PTAG(Formulation)), 0> type;
};

/*!
 * \brief Set the property for the material law by retrieving it from
 *        the spatial parameters.
 */
SET_PROP(BoxTwoP, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;

public:
    typedef typename SpatialParameters::MaterialLaw type;
};

/*!
 * \brief Set the property for the material parameters by extracting
 *        it from the material law.
 */
SET_PROP(BoxTwoP, MaterialLawParams)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLaw)) MaterialLaw;

public:
    typedef typename MaterialLaw::Params type;
};

SET_TYPE_PROP(BoxTwoP, FluidSystem, FluidSystem2P<TypeTag>);

SET_TYPE_PROP(BoxTwoP, FluidState, TwoPFluidState<TypeTag>);

// enable jacobian matrix recycling by default
SET_BOOL_PROP(BoxTwoP, EnableJacobianRecycling, true);
// enable partial reassembling by default
SET_BOOL_PROP(BoxTwoP, EnablePartialReassemble, true);
// enable time-step ramp up by default
SET_BOOL_PROP(BoxTwoP, EnableTimeStepRampUp, false);

}
//

}

#endif
