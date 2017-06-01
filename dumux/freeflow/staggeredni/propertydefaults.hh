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
 * \ingroup OnePModel
 * \file
 *
 * \brief Defines the properties required for the one-phase fully implicit model.
 */
#ifndef DUMUX_NAVIER_STOKES_NI_PROPERTY_DEFAULTS_HH
#define DUMUX_NAVIER_STOKES_NI_PROPERTY_DEFAULTS_HH

#include "properties.hh"
#include "model.hh"
#include "../staggered/volumevariables.hh"
#include "indices.hh"
#include "../staggered/propertydefaults.hh" //TODO: why do we need this include?

namespace Dumux
{

// \{

///////////////////////////////////////////////////////////////////////////
// default property values for the non-isothermal single phase model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

SET_PROP(StaggeredNonIsothermal, NumEqCellCenter)
{
private:
    static constexpr auto isothermalNumEqCellCenter = GET_PROP_VALUE(TypeTag, IsothermalNumEqCellCenter);
public:
    static constexpr auto value = isothermalNumEqCellCenter + 1;
};

// SET_INT_PROP(StaggeredNonIsothermal, NumEqCellCenter, 2);
//
// //! the VolumeVariables property
SET_TYPE_PROP(StaggeredNonIsothermal, VolumeVariables, NavierStokesVolumeVariables<TypeTag>);
SET_TYPE_PROP(StaggeredNonIsothermal, Model, StaggeredNonIsothermalModel<TypeTag>);
SET_TYPE_PROP(StaggeredNonIsothermal, Indices, StaggeredNonIsothermalIndices<TypeTag>);
//
SET_BOOL_PROP(StaggeredNonIsothermal, EnableEnergyBalanceStokes, true);
//
SET_BOOL_PROP(StaggeredNonIsothermal, UseMoles, true);
//
SET_TYPE_PROP(StaggeredNonIsothermal, HeatConductionType, FouriersLaw<TypeTag>);
//
SET_INT_PROP(StaggeredNonIsothermal, PhaseIdx, 0); //!< Defines the phaseIdx


} // end namespace Properties

} // end namespace Dumux

#endif
