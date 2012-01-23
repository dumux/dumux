// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010 by Markus Wolff                                      *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
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
 * \ingroup OnePhase
 * \ingroup Properties
 * \file
 *
 * \brief Defines the properties required for the single phase sequential model.
 */

#ifndef DUMUX_1PPROPERTIES_HH
#define DUMUX_1PPROPERTIES_HH

//Dumux-includes
#include <dumux/decoupled/common/decoupledproperties.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include "1pindices.hh"

namespace Dumux
{

////////////////////////////////
// forward declarations
////////////////////////////////

////////////////////////////////
// properties
////////////////////////////////
namespace Properties
{
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the single-phase problem
NEW_TYPE_TAG(DecoupledOneP, INHERITS_FROM(DecoupledModel));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG( SpatialParameters )
; //!< The type of the spatial parameters object
NEW_PROP_TAG( EnableGravity)
; //!< Returns whether gravity is considered in the problem
NEW_PROP_TAG( Fluid )
; //!< The fluid for one-phase models
}
}

#include <dumux/decoupled/common/variableclass.hh>
#include <dumux/decoupled/1p/cellData1p.hh>

namespace Dumux
{
namespace Properties
{
//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

SET_INT_PROP(DecoupledOneP, NumEq, 1);

SET_INT_PROP(DecoupledOneP, NumPhases, 1)//!< Single phase system
;
SET_INT_PROP(DecoupledOneP, NumComponents, 1); //!< Each phase consists of 1 pure component

SET_TYPE_PROP(DecoupledOneP, Indices, DecoupledOnePCommonIndices);

SET_TYPE_PROP(DecoupledOneP, Variables, VariableClass<TypeTag>);

SET_TYPE_PROP(DecoupledOneP, CellData, CellData1P<TypeTag>);
}
}
#endif
