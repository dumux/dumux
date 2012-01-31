// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Markus Wolff                                      *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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
#ifndef DUMUX_IMPES2PADAPTIVE_PROPERTIES_HH
#define DUMUX_IMPES2PADAPTIVE_PROPERTIES_HH

#include <dumux/decoupled/common/impetproperties.hh>
#include <dumux/decoupled/2p/2pproperties.hh>

/*!
 * \ingroup IMPET
 * \ingroup Properties
 */
/*!
 * \file
 * \brief Base file for properties related to sequential IMPET algorithms
 */
namespace Dumux
{

namespace Properties
{
/*!
 *
 * \brief General properties for sequential IMPET algorithms
 *
 * This class holds properties necessary for the sequential IMPET solution.
 */

//////////////////////////////////////////////////////////////////
// Type tags tags
//////////////////////////////////////////////////////////////////

//! The type tag for models based on the diffusion-scheme
NEW_TYPE_TAG(IMPESTwoPAdaptive, INHERITS_FROM(IMPET, DecoupledTwoP));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
}
}

#include <dumux/decoupled/common/variableclass_adaptive.hh>
#include <dumux/decoupled/2p/cellData2p_adaptive.hh>
#include "gridadaptionindicator2p.hh"

namespace Dumux
{
namespace Properties
{
SET_BOOL_PROP(IMPESTwoPAdaptive, AdaptiveGrid, true);
SET_TYPE_PROP(IMPESTwoPAdaptive, Variables, Dumux::VariableClassAdaptive<TypeTag>);
SET_TYPE_PROP(IMPESTwoPAdaptive, CellData, Dumux::CellData2PAdaptive<TypeTag>);
SET_TYPE_PROP(IMPESTwoPAdaptive, AdaptionIndicator, GridAdaptionIndicator2P<TypeTag>);
}
}

#endif
