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
 * \ingroup NIModel
 * \file
 *
 * \brief Defines default values for most properties required by the
 *        implicit non-isothermal models.
 */
#ifndef DUMUX_MIMETIC_ENERGY_PROPERTY_DEFAULTS_HH
#define DUMUX_MIMETIC_ENERGY_PROPERTY_DEFAULTS_HH

#include "indices.hh"

namespace Dumux
{

namespace Properties
{

//////////////////////////////////////////////////////////////////
// Property values
//////////////////////////////////////////////////////////////////

SET_BOOL_PROP(MimeticNonIsothermal, EnableEnergyBalance, true);

//! add the energy balance
//SET_INT_PROP(MimeticNonIsothermal, NumEq, GET_PROP_VALUE(TypeTag, IsothermalNumEq) + 1);

SET_INT_PROP(MimeticNonIsothermal, NumEqCellCenter, GET_PROP_VALUE(TypeTag, IsothermalNumEqCellCenter) + 1);
SET_INT_PROP(MimeticNonIsothermal, NumEqFace,  GET_PROP_VALUE(TypeTag, IsothermalNumEqFace) + 1);

//! indices for non-isothermal models
SET_TYPE_PROP(MimeticNonIsothermal, Indices, MimeticEnergyIndices<TypeTag, 0>);

} // end namespace Properties

} // end namespace Dumux

#endif
