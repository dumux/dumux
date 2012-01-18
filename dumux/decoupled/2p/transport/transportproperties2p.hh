// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Markus Wolff                                      *
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
#ifndef DUMUX_TRANSPORT_PROPERTIES_2P_HH
#define DUMUX_TRANSPORT_PROPERTIES_2P_HH

#include <dumux/decoupled/common/transportproperties.hh>

/*!
 * \ingroup Saturation2p
 * \ingroup Properties
 */
/*!
 * \file
 * \brief Specifies the properties for immiscible 2p transport
 */
namespace Dumux
{

template<class TypeTag>
class DiffusivePart;

template<class TypeTag>
class ConvectivePart;

template<class TypeTag>
class EvalCflFluxDefault;

template<class TypeTag>
class FVVelocityDefault;

namespace Properties
{
// \{

//////////////////////////////////////////////////////////////////
// Type tags tags
//////////////////////////////////////////////////////////////////

//! The type tag for models based on the diffusion-scheme
NEW_TYPE_TAG(TransportTwoP, INHERITS_FROM(Transport));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG( CapillaryFlux );         //!< The type of the diffusive part in a transport equation
NEW_PROP_TAG( GravityFlux );        //!< The type of a convective part in a transport equation
NEW_PROP_TAG( CalculateVelocityInTransport );
NEW_PROP_TAG( PrecomputedConstRels );

SET_TYPE_PROP(TransportTwoP, CapillaryFlux, DiffusivePart<TypeTag>);
SET_TYPE_PROP(TransportTwoP, GravityFlux, ConvectivePart<TypeTag>);
SET_TYPE_PROP(TransportTwoP, EvalCflFluxFunction, EvalCflFluxDefault<TypeTag>);
SET_TYPE_PROP(TransportTwoP, Velocity, FVVelocityDefault<TypeTag>);
SET_BOOL_PROP( TransportTwoP, CalculateVelocityInTransport, true);
SET_BOOL_PROP( TransportTwoP, PrecomputedConstRels, true);
}
}

#endif
