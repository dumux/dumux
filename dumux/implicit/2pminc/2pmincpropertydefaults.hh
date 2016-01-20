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
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/
/*!
 * \ingroup TwoPTwoCModel
 */
/*!
 * \file
 *
 * \brief Defines additional default values for the properties
 *        required for the coupling of the 2p2cni box model
 */
#ifndef DUMUX_2PMIN_PROPERTY_DEFAULTS_HH
#define DUMUX_2PMIN_PROPERTY_DEFAULTS_HH

#include <dumux/porousmediumflow/2p/implicit/propertydefaults.hh>

#include "2pmincproperties.hh"
#include "2pmincmodel.hh"
#include "2pminclocalresidual.hh"
#include "2pmincvolumevariables.hh"
#include "2pmincfluxvariables.hh"
#include "2pmincindices.hh"

namespace Dumux
{

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Property values
//////////////////////////////////////////////////////////////////
//! Use the 2pminc local jacobian operator for the 2pminc model
SET_TYPE_PROP(TwoPMinc,
              LocalResidual,
              TwoPMincLocalResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(TwoPMinc, Model, TwoPMincModel<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(TwoPMinc, VolumeVariables, TwoPMincVolumeVariables<TypeTag>);
//! the FluxVariables property
SET_TYPE_PROP(TwoPMinc, FluxVariables, TwoPMincFluxVariables<TypeTag>);

//! The Indices property
SET_TYPE_PROP(TwoPMinc,
              Indices,
              TwoPMincIndices<TypeTag, GET_PROP_VALUE(TypeTag, Formulation), 0>);


//! Set the number of continua to 2
SET_INT_PROP(TwoPMinc, NumContinua, 2);

//! Set the default inteaction type to constant volume fractions
SET_INT_PROP(TwoPMinc, ProblemInteractingContinuaType, 0);
}
}
#endif
