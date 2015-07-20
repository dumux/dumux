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
 * \brief Defines the properties required for the ZeroEq model.
 */

#ifndef DUMUX_ZEROEQ_PROPERTY_DEFAULTS_HH
#define DUMUX_ZEROEQ_PROPERTY_DEFAULTS_HH

#include "zeroeqfluxvariables.hh"
#include "zeroeqmodel.hh"

namespace Dumux
{

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

// INHERITED
//! Set the Model property
SET_TYPE_PROP(BoxZeroEq, Model, ZeroEqModel<TypeTag>);

//! Set the Flux Variables property
SET_TYPE_PROP(BoxZeroEq, FluxVariables, ZeroEqFluxVariables<TypeTag>);

//! Set the indices used by the ZeroEq model
SET_TYPE_PROP(BoxZeroEq, Indices, ZeroEqCommonIndices<TypeTag>);

//! Set calculation to Navier Stokes
SET_BOOL_PROP(BoxZeroEq, EnableNavierStokes, true);

// NEW PROPERTIES
//! Set the Karman constant \f$[-]\f$
SET_SCALAR_PROP(BoxZeroEq, KarmanConstant, 0.41);

//! Set BBoxMin of wall normal direction as a wall
SET_BOOL_PROP(BoxZeroEq, BBoxMinIsWall, true);

//! Set BBoxMax of wall normal direction as a wall
SET_BOOL_PROP(BoxZeroEq, BBoxMaxIsWall, true);

//! Set main flow direction
SET_INT_PROP(BoxZeroEq, ZeroEqFlowNormal, 0);

//! Set wall normal direction
SET_INT_PROP(BoxZeroEq, ZeroEqWallNormal, 1);

//! Set a zero sand grain roughness \f$[m]\f$ at BBoxMin
SET_SCALAR_PROP(BoxZeroEq, ZeroEqBBoxMinSandGrainRoughness, 0.0);

//! Set a zero sand grain roughness \f$[m]\f$ at BBoxMax
SET_SCALAR_PROP(BoxZeroEq, ZeroEqBBoxMaxSandGrainRoughness, 0.0);

//! Set the eddy viscosity model
SET_INT_PROP(BoxZeroEq, ZeroEqEddyViscosityModel, 1);

//! Set the number of wall intervals in which more complex turbulence models are evaluated
SET_INT_PROP(BoxZeroEq, NumberOfIntervals, 1000);

//! Set the BaseStokesFluxVariables to StokesFluxVariables
SET_TYPE_PROP(BoxZeroEq, BaseStokesFluxVariables, StokesFluxVariables<TypeTag>);

//! Set the BaseStokesModel to StokesModel
SET_TYPE_PROP(BoxZeroEq, BaseStokesModel, StokesModel<TypeTag>);
}
} // namespace Dumux
#endif // DUMUX_ZEROEQ_PROPERTY_DEFAULTS_HH
