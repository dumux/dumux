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
 * \ingroup IMPES
 * \ingroup IMPETProperties
 */
/*!
 * \file
 * \brief Properties for adaptive implementations of the sequential IMPES algorithms
 */
namespace Dumux
{

namespace Properties
{
/*!
 *
 * \brief General properties for adaptive implementations of the sequential IMPES algorithms
 */

//////////////////////////////////////////////////////////////////
// Type tags tags
//////////////////////////////////////////////////////////////////

//!  TypeTag for grid-adaptive two-phase IMPES scheme
NEW_TYPE_TAG(IMPESTwoPAdaptive, INHERITS_FROM(IMPET, DecoupledTwoP));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
}
}

#include <dumux/decoupled/common/variableclassadaptive.hh>
#include <dumux/decoupled/2p/cellData2padaptive.hh>
#include "gridadaptionindicator2p.hh"
#include <dumux/decoupled/2p/impes/impesproblem2p.hh>
#include <dumux/decoupled/common/gridadaptinitializationindicator.hh>

namespace Dumux
{
namespace Properties
{
//! Enable adaptive grid
SET_BOOL_PROP(IMPESTwoPAdaptive, AdaptiveGrid, true);
//! Set variable class for adaptive impet schemes
SET_TYPE_PROP(IMPESTwoPAdaptive, Variables, Dumux::VariableClassAdaptive<TypeTag>);
//! Set cell data class for adaptive two-phase IMPES schemes
SET_TYPE_PROP(IMPESTwoPAdaptive, CellData, Dumux::CellData2PAdaptive<TypeTag>);
//! Set the standard indicator class of two-phase models for adaption or coarsening
SET_TYPE_PROP(IMPESTwoPAdaptive, AdaptionIndicator, GridAdaptionIndicator2P<TypeTag>);
//!Set default class for adaptation initialization indicator
SET_TYPE_PROP(IMPESTwoPAdaptive,  AdaptionInitializationIndicator, GridAdaptInitializationIndicator<TypeTag>);
}
}

#endif
