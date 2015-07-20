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
 * \ingroup BoxZeroEqModel
 *
 * \file
 *
 * \brief Defines the supplementary properties required for the
 *        ZeroEq model.
 *
 */

#ifndef DUMUX_ZEROEQ_PROPERTIES_HH
#define DUMUX_ZEROEQ_PROPERTIES_HH

#include <dumux/freeflow/stokes/stokesproperties.hh>

namespace Dumux
{

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the ZeroEq Problem
NEW_TYPE_TAG(BoxZeroEq, INHERITS_FROM(BoxStokes));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(KarmanConstant); //!< The Karman constant
NEW_PROP_TAG(BBoxMinIsWall); //!< Sets BBoxMin as a wall
NEW_PROP_TAG(BBoxMaxIsWall); //!< Sets BBoxMax as a wall
NEW_PROP_TAG(ZeroEqFlowNormal); //!< Indicates the main flow direction
NEW_PROP_TAG(ZeroEqWallNormal); //!< Indicates the wall normal direction
NEW_PROP_TAG(ZeroEqBBoxMinSandGrainRoughness); //!< Sets a sand grain roughness at BBoxMin
NEW_PROP_TAG(ZeroEqBBoxMaxSandGrainRoughness); //!< Sets a sand grain roughness at BBoxMax
NEW_PROP_TAG(ZeroEqEddyViscosityModel); //!< Returns used eddy viscosity model
NEW_PROP_TAG(NumberOfIntervals); //!< Returns number of wall intervals
NEW_PROP_TAG(BaseStokesModel); //!< Returns the base implementation of the Stokes model
NEW_PROP_TAG(BaseStokesVolumeVariables); //!< Returns the base implementation of the Stokes volume variables
NEW_PROP_TAG(BaseStokesFluxVariables); //!< Returns the base implementation of the Stokes flux variables
}

}

#endif // DUMUX_ZEROEQ_PROPERTIES_HH
