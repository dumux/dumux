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
 * \file
 * \ingroup ZeroEqModel
 *
 * \brief A single-phase, isothermal Reynolds-Averaged Navier-Stokes 0-Eq. model
 *
 * \copydoc RANSModel
 *
 * These models calculate the eddy viscosity without solving additional PDEs,
 * only based on the wall distance and the velocity gradient.
 *
 * \copydoc Dumux::EddyViscosityModels
 */

#ifndef DUMUX_ZEROEQ_MODEL_HH
#define DUMUX_ZEROEQ_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/freeflow/properties.hh>
#include <dumux/freeflow/rans/model.hh>

#include "volumevariables.hh"

namespace Dumux {
namespace Properties {

///////////////////////////////////////////////////////////////////////////
// default property values for the isothermal RANS 0-Eq. model
///////////////////////////////////////////////////////////////////////////

//! The type tag for the single-phase, isothermal Reynolds-Averaged Navier-Stokes 0-Eq. model
NEW_TYPE_TAG(ZeroEq, INHERITS_FROM(RANS));

//! The volume variables
SET_TYPE_PROP(ZeroEq, VolumeVariables, ZeroEqVolumeVariables<TypeTag>);

//////////////////////////////////////////////////////////////////
// default property values for the non-isothermal RANS 0-Eq. model
//////////////////////////////////////////////////////////////////

//! The type tag for the single-phase, isothermal Reynolds-Averaged Navier-Stokes model
NEW_TYPE_TAG(ZeroEqNI, INHERITS_FROM(RANSNI));

//! The volume variables
SET_TYPE_PROP(ZeroEqNI, VolumeVariables, ZeroEqVolumeVariables<TypeTag>);

// \}
}

} // end namespace

#endif // DUMUX_ZEROEQ_MODEL_HH
