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
#ifndef DUMUX_FV_TRANSPORT_PROPERTIES_2P_HH
#define DUMUX_FV_TRANSPORT_PROPERTIES_2P_HH

#include <dumux/decoupled/2p/transport/transportproperties2p.hh>

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
//////////////////////////////////////////////////////////////////
// Forward Declarations
//////////////////////////////////////////////////////////////////
template <class TypeTag>
class FVSaturation2P;

namespace Properties
{
// \{

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for models based on the diffusion-scheme
NEW_TYPE_TAG(FVTransportTwoP, INHERITS_FROM(TransportTwoP));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG( CapillaryFlux );         //!< The type of the diffusive part in a transport equation
NEW_PROP_TAG( GravityFlux );        //!< The type of a convective part in a transport equation
NEW_PROP_TAG( PrecomputedConstRels );
}
}

#include "evalcflflux_default.hh"
#include "dumux/decoupled/2p/transport/fv/diffusivepart.hh"
#include "dumux/decoupled/2p/transport/fv/convectivepart.hh"
#include <dumux/decoupled/2p/transport/fv/fvsaturation2p.hh>

namespace Dumux
{
namespace Properties
{
SET_TYPE_PROP(FVTransportTwoP, EvalCflFluxFunction, EvalCflFluxDefault<TypeTag>);
SET_TYPE_PROP(FVTransportTwoP, CapillaryFlux, DiffusivePart<TypeTag>);
SET_TYPE_PROP(FVTransportTwoP, GravityFlux, ConvectivePart<TypeTag>);
SET_BOOL_PROP( FVTransportTwoP, PrecomputedConstRels, true);

// Set the model properties
SET_PROP(FVTransportTwoP, TransportModel)
{
    typedef Dumux::FVSaturation2P<TypeTag> type;
};
}
}

#endif
