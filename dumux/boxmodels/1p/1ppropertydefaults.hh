// $Id$
/*****************************************************************************
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
 * \file
 *
 * \brief Defines the properties required for the one-phase BOX model.
 */
#ifndef DUMUX_1P_PROPERTY_DEFAULTS_HH
#define DUMUX_1P_PROPERTY_DEFAULTS_HH

#include <dumux/boxmodels/common/boxproperties.hh>

#include "1pmodel.hh"
#include "1plocalresidual.hh"
#include "1pvolumevariables.hh"
#include "1pfluxvariables.hh"
#include "1pindices.hh"

namespace Dumux
{
/*!
 * \addtogroup OnePBoxModel
 */
// \{

///////////////////////////////////////////////////////////////////////////
// default property values for the isothermal single phase model
///////////////////////////////////////////////////////////////////////////
namespace Properties {
SET_INT_PROP(BoxOneP, NumEq, 1); //!< set the number of equations to 1
SET_INT_PROP(BoxOneP, NumPhases, 1); //!< The number of phases in the 1p model is 1

//! The local residual function
SET_TYPE_PROP(BoxOneP,
              LocalResidual,
              OnePLocalResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(BoxOneP, Model, OnePBoxModel<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(BoxOneP, VolumeVariables, OnePVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(BoxOneP, FluxVariables, OnePFluxVariables<TypeTag>);

//! The indices required by the isothermal single-phase model
SET_TYPE_PROP(BoxOneP, OnePIndices, OnePIndices);

//! The weight of the upwind control volume when calculating
//! fluxes. Use central differences by default.
SET_SCALAR_PROP(BoxOneP, UpwindWeight, 0.5);

// \}
};

} // end namepace

#endif
