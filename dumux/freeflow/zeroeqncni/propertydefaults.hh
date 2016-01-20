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
 * \ingroup BoxZeroEqcniModel
 * \file
 *
 * \brief Defines default properties for the non-isothermal compositional ZeroEq box model.
 */

#ifndef DUMUX_ZEROEQNCNI_PROPERTY_DEFAULTS_HH
#define DUMUX_ZEROEQNCNI_PROPERTY_DEFAULTS_HH

#include <dumux/freeflow/stokesncni/fluxvariables.hh>
#include <dumux/freeflow/stokesncni/model.hh>
#include <dumux/freeflow/stokesncni/volumevariables.hh>
#include "fluxvariables.hh"
#include "model.hh"

namespace Dumux
{

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

// INHERITED
//! Set the Model property
SET_TYPE_PROP(BoxZeroEqncni, Model, ZeroEqncniModel<TypeTag>);

//! Set the FluxVariables property
SET_TYPE_PROP(BoxZeroEqncni, FluxVariables, ZeroEqncniFluxVariables<TypeTag>);

//! Set the indices used by the ZeroEqncni model
SET_TYPE_PROP(BoxZeroEqncni, Indices, ZeroEqncniCommonIndices<TypeTag>);

//! Set calculation to Navier-Stokes
SET_BOOL_PROP(BoxZeroEqncni, EnableNavierStokes, true);

//! Don't use mole fractions
SET_BOOL_PROP(BoxZeroEqncni, UseMoles, false);

// NEW PROPERTIES
//! Set the Eddy Conductivity Model
SET_INT_PROP(BoxZeroEqncni, ZeroEqEddyConductivityModel, 1);

//! Set the turbulent Prandtl number \f$[-]\f$ (Reynolds analogy)
SET_SCALAR_PROP(BoxZeroEqncni, ZeroEqTurbulentPrandtlNumber, 1.0);

//! Set the BaseStokesModel to StokesncniModel
SET_TYPE_PROP(BoxZeroEqncni, BaseStokesModel, StokesncniModel<TypeTag>);

//! Set the BaseStokesFluxVariables to StokesncniFluxVariables
SET_TYPE_PROP(BoxZeroEqncni, BaseStokesFluxVariables, StokesncniFluxVariables<TypeTag>);
}
}

#endif // DUMUX_ZEROEQNCNI_PROPERTY_DEFAULTS_HH
