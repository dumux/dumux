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
#ifndef DUMUX_TRANSPORT_PROPERTIES_HH
#define DUMUX_TRANSPORT_PROPERTIES_HH

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

namespace Properties
{
// \{

//////////////////////////////////////////////////////////////////
// Type tags tags
//////////////////////////////////////////////////////////////////

//! The type tag for models based on the diffusion-scheme
NEW_TYPE_TAG(Transport);

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG( DiffusivePart );         //!< The type of the diffusive part in a transport equation
NEW_PROP_TAG( ConvectivePart );        //!< The type of a convective part in a transport equation
NEW_PROP_TAG( EvalCflFluxFunction ); //!< Type of the evaluation of the CFL-condition
NEW_PROP_TAG( CFLFactor );

SET_TYPE_PROP(Transport, DiffusivePart, DiffusivePart<TypeTag>);
SET_TYPE_PROP(Transport, ConvectivePart, ConvectivePart<TypeTag>);
SET_TYPE_PROP(Transport, EvalCflFluxFunction, EvalCflFluxDefault<TypeTag>);
SET_SCALAR_PROP(Transport, CFLFactor, 1.0);
}
}

#endif
