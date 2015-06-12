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
 * \ingroup BoxZeroEqncModel
 *
 * \file
 *
 * \brief Defines default properties for the compositional ZeroEq box model.
 */

#ifndef DUMUX_ZEROEQNC_PROPERTY_DEFAULTS_HH
#define DUMUX_ZEROEQNC_PROPERTY_DEFAULTS_HH

#include <dumux/freeflow/stokesnc/stokesncfluxvariables.hh>
#include <dumux/freeflow/stokesnc/stokesncmodel.hh>
#include <dumux/freeflow/stokesnc/stokesncvolumevariables.hh>
#include "zeroeqncfluxvariables.hh"
#include "zeroeqncmodel.hh"

namespace Dumux
{

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

// INHERITED
//! Set the Model property
SET_TYPE_PROP(BoxZeroEqnc, Model, ZeroEqncModel<TypeTag>);

//! Set the FluxVariables property
SET_TYPE_PROP(BoxZeroEqnc, FluxVariables, ZeroEqncFluxVariables<TypeTag>);

//! Set the indices used by the ZeroEqnc model
SET_TYPE_PROP(BoxZeroEqnc, Indices, ZeroEqncCommonIndices<TypeTag>);

//! Set calculation to Navier-Stokes
SET_BOOL_PROP(BoxZeroEqnc, EnableNavierStokes, true);

//! Don't use mole fractions
SET_BOOL_PROP(BoxZeroEqnc, UseMoles, false);

// NEW PROPERTIES
//! Set the Eddy Diffusivity Model
SET_INT_PROP(BoxZeroEqnc, ZeroEqEddyDiffusivityModel, 1);

//! Set the turbulent Schmidt number \f$[-]\f$ (Reynolds analogy)
SET_SCALAR_PROP(BoxZeroEqnc, ZeroEqTurbulentSchmidtNumber, 1.0);

//! Set the BaseStokesModel to StokesncModel
SET_TYPE_PROP(BoxZeroEqnc, BaseStokesModel, StokesncModel<TypeTag>);

//! Set the BaseStokesFluxVariables to StokesncniFluxVariables
SET_TYPE_PROP(BoxZeroEqnc, BaseStokesFluxVariables, StokesncFluxVariables<TypeTag>);

}
}

#endif // DUMUX_ZEROEQNC_PROPERTY_DEFAULTS_HH
