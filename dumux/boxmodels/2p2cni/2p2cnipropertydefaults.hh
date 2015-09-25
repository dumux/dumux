// $Id$
/*****************************************************************************
 *   Copyright (C) 2008-2010 by Andreas Lauser                               *
 *   Copyright (C) 2008-2009 by Melanie Darcis                               *
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2009 by Bernd Flemisch                               *
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
 * \ingroup TwoPTwoCNIModel
 */
/*!
 * \file
 *
 * \brief Defines default values for most properties required by the 2p2cni
 *        box model.
 */
#ifndef DUMUX_2P2CNI_PROPERTY_DEFAULTS_HH
#define DUMUX_2P2CNI_PROPERTY_DEFAULTS_HH

#include <dumux/boxmodels/2p2c/2p2cpropertydefaults.hh>

#include "2p2cnimodel.hh"
#include "2p2cniproblem.hh"
#include "2p2cniindices.hh"
#include "2p2cnilocalresidual.hh"
#include "2p2cnivolumevariables.hh"
#include "2p2cnifluxvariables.hh"

namespace Dumux
{

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Property values
//////////////////////////////////////////////////////////////////

SET_INT_PROP(BoxTwoPTwoCNI, NumEq, 3); //!< set the number of equations to 3

//! Use the 2p2cni local jacobian operator for the 2p2cni model
SET_TYPE_PROP(BoxTwoPTwoCNI,
              LocalResidual,
              TwoPTwoCNILocalResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(BoxTwoPTwoCNI, Model, TwoPTwoCNIModel<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(BoxTwoPTwoCNI, VolumeVariables, TwoPTwoCNIVolumeVariables<TypeTag>);


//! the FluxVariables property
SET_TYPE_PROP(BoxTwoPTwoCNI, FluxVariables, TwoPTwoCNIFluxVariables<TypeTag>);

//! The indices required by the non-isothermal 2p2c model
SET_PROP(BoxTwoPTwoCNI, TwoPTwoCIndices)
{ typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPTwoCNIIndices)) type; };

SET_PROP(BoxTwoPTwoCNI, TwoPTwoCNIIndices)
{ private:
    enum { formulation = GET_PROP_VALUE(TypeTag, PTAG(Formulation)) };
public:
    typedef TwoPTwoCNIIndices<TypeTag, formulation, 0> type;
};

}

}
#endif
