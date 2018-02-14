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
 * \ingroup SweModel
 *
 * \brief A two-dimesnional shallow water equations model
 *
 *
 *
 * So far, only the cell centerd spatial discretization is available.
 */

#ifndef DUMUX_SWE_MODEL_HH
#define DUMUX_SWE_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/shallowwater/properties.hh>

#include "localresidual.hh"
#include "volumevariables.hh"
#include "fluxvariables.hh"
#include "fluxvariablescache.hh"
#include "indices.hh"
#include "vtkoutputfields.hh"

#include <dumux/discretization/methods.hh>

namespace Dumux
{

///////////////////////////////////////////////////////////////////////////
// properties for the shallow water model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the Swe inherits from shallow water
NEW_TYPE_TAG(Swe, INHERITS_FROM(ModelProperties));

///////////////////////////////////////////////////////////////////////////
// default property values for shallow water model
///////////////////////////////////////////////////////////////////////////
SET_BOOL_PROP(Swe, EnableNumericalFlux, true); //!< Enable numerical flux (e.g. Riemann solver)
SET_BOOL_PROP(Swe, EnableTurbulenceModel, false); //!< The turbulence model for the momentum equations


/*!
* \brief The two-dimensional shallow water equations have allways 3 equations.
*/
SET_INT_PROP(Swe, NumEq, 3);
SET_INT_PROP(Swe, NumEqVector, 3);

//! The local residual
SET_TYPE_PROP(Swe, LocalResidual, SweResidual<TypeTag>);

//! The volume variables
SET_TYPE_PROP(Swe, VolumeVariables, SweVolumeVariables<TypeTag>);

//! The flux variables
SET_TYPE_PROP(Swe, FluxVariables, SweFluxVariables<TypeTag>);

//! The flux variables cache class
SET_TYPE_PROP(Swe, FluxVariablesCache, SweFluxVariablesCache<TypeTag>);

//! The indices required for the SWEs
SET_TYPE_PROP(Swe, Indices, SweIndices<TypeTag>);

//! The specific vtk output fields
SET_TYPE_PROP(Swe, VtkOutputFields, SweVtkOutputFields<TypeTag>);

}

} // end namespace

#endif // DUMUX_SWE_MODEL_HH
