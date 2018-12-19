// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \file
 * \ingroup SequentialTwoPModel
 * \brief Properties for adaptive implementations of the sequential IMPES algorithms
 */
#ifndef DUMUX_IMPES2PADAPTIVE_PROPERTIES_HH
#define DUMUX_IMPES2PADAPTIVE_PROPERTIES_HH

#include <dumux/porousmediumflow/sequential/impetproperties.hh>
#include <dumux/porousmediumflow/2p/sequential/properties.hh>

namespace Dumux {
namespace Properties {
/*!
 * \ingroup SequentialTwoPModel
 * \brief General properties for adaptive implementations of the sequential IMPES algorithms
 */

//////////////////////////////////////////////////////////////////
// Type tags tags
//////////////////////////////////////////////////////////////////

//!  TypeTag for grid-adaptive two-phase IMPES scheme
NEW_TYPE_TAG(IMPESTwoPAdaptive, INHERITS_FROM(IMPET, SequentialTwoP));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
} // end namespace Properties
} // end namespace Dumux

#include <dumux/porousmediumflow/sequential/variableclassadaptive.hh>
#include <dumux/porousmediumflow/2p/sequential/celldataadaptive.hh>
#include "gridadaptionindicator.hh"
#include <dumux/porousmediumflow/2p/sequential/impes/problem.hh>
#include <dumux/porousmediumflow/sequential/gridadaptinitializationindicator.hh>

namespace Dumux {
namespace Properties {
//! Enable adaptive grid
SET_BOOL_PROP(IMPESTwoPAdaptive, AdaptiveGrid, true);
//! Set variable class for adaptive impet schemes
SET_TYPE_PROP(IMPESTwoPAdaptive, Variables, VariableClassAdaptive<TypeTag>);
//! Set cell data class for adaptive two-phase IMPES schemes
SET_TYPE_PROP(IMPESTwoPAdaptive, CellData, CellData2PAdaptive<TypeTag>);
//! Set the standard indicator class of two-phase models for adaption or coarsening
SET_TYPE_PROP(IMPESTwoPAdaptive, AdaptionIndicator, GridAdaptionIndicator2P<TypeTag>);
//! Set default class for adaptation initialization indicator
SET_TYPE_PROP(IMPESTwoPAdaptive,  AdaptionInitializationIndicator, GridAdaptInitializationIndicator<TypeTag>);
} // end namespace Properties
} // end namespace Dumux

#endif
