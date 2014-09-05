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
 * \ingroup TwoPTwoCNIModel
 * \file
 *
 * \brief Defines default values for most properties required by the
 *        non-isothermal two-phase two-component fully implicit model.
 */
#ifndef DUMUX_2P2CNI_PROPERTY_DEFAULTS_HH
#define DUMUX_2P2CNI_PROPERTY_DEFAULTS_HH

#include "2p2cniproperties.hh"
#include <dumux/implicit/2p2c/2p2cpropertydefaults.hh>

#include "2p2cnimodel.hh"
#include "2p2cniindices.hh"
#include "2p2cnilocalresidual.hh"
#include "2p2cnivolumevariables.hh"
#include "2p2cnifluxvariables.hh"

#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivitysomerton.hh>

namespace Dumux
{

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Property values
//////////////////////////////////////////////////////////////////
//! Set the number of equations to 3
SET_INT_PROP(TwoPTwoCNI, NumEq, 3); //!< set the number of equations to 3

//! Use the 2p2cni local jacobian operator for the 2p2cni model
SET_TYPE_PROP(TwoPTwoCNI,
              LocalResidual,
              TwoPTwoCNILocalResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(TwoPTwoCNI, Model, TwoPTwoCNIModel<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(TwoPTwoCNI, VolumeVariables, TwoPTwoCNIVolumeVariables<TypeTag>);


//! the FluxVariables property
SET_TYPE_PROP(TwoPTwoCNI, FluxVariables, TwoPTwoCNIFluxVariables<TypeTag>);

//! The indices required by the non-isothermal 2p2c model
SET_PROP(TwoPTwoCNI, Indices) 
{ private:
    enum { formulation = GET_PROP_VALUE(TypeTag, Formulation) };
 public:
    typedef TwoPTwoCNIIndices<TypeTag, formulation, 0> type;
};

//! Somerton is used as default model to compute the effective thermal heat conductivity
SET_PROP(TwoPTwoCNI, ThermalConductivityModel)
{ private :
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
  public:
    typedef ThermalConductivitySomerton<Scalar> type;
};

}

}
#endif
