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
#ifndef DUMUX_NAVIER_STOKES_NI_PROPERTIES_HH
#define DUMUX_NAVIER_STOKES_NI_PROPERTIES_HH

#include <dumux/common/properties.hh>
#include "fluxvariables.hh"
#include "indices.hh"
#include "localresidual.hh"
#include "vtkoutputfields.hh"
#include <dumux/discretization/fourierslaw.hh>

namespace Dumux
{

// \{

namespace Properties {

//! The type tags for the non-isothermal Navier Stokes problems
NEW_TYPE_TAG(NavierStokesNonIsothermal);

NEW_PROP_TAG(IsothermalNumEqCellCenter);
NEW_PROP_TAG(IsothermalNumEqFace);

///////////////////////////////////////////////////////////////////////////
// default property values for the non-isothermal single phase model
///////////////////////////////////////////////////////////////////////////

SET_PROP(NavierStokesNonIsothermal, NumEqCellCenter)
{
private:
    static constexpr auto isothermalNumEqCellCenter = GET_PROP_VALUE(TypeTag, IsothermalNumEqCellCenter);
public:
    static constexpr auto value = isothermalNumEqCellCenter + 1;
};

SET_PROP(NavierStokesNonIsothermal, NumEq)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    static constexpr auto dim = GridView::dimension;
    static constexpr auto numComponents = GET_PROP_VALUE(TypeTag, NumComponents);
public:
    static constexpr int value = dim + numComponents + 1;
};

SET_TYPE_PROP(NavierStokesNonIsothermal, Indices, NavierStokesNonIsothermalIndices<TypeTag>);

SET_TYPE_PROP(NavierStokesNonIsothermal, VtkOutputFields, FreeFlowEnergyVtkOutputFields<TypeTag>);

SET_BOOL_PROP(NavierStokesNonIsothermal, EnableEnergyBalance, true);

SET_TYPE_PROP(NavierStokesNonIsothermal, HeatConductionType, FouriersLaw<TypeTag>);

} // end namespace Properties

} // end namespace Dumux

#endif
