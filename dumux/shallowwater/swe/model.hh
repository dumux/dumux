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
#include <dumux/common/properties/model.hh>
#include <dumux/shallowwater/properties.hh>
#include <dumux/shallowwater/advectiveflux.hh>

#include "localresidual.hh"
#include "volumevariables.hh"
#include "fluxvariables.hh"
#include "fluxvariablescache.hh"
#include "indices.hh"
#include "vtkoutputfields.hh"

#include <dumux/discretization/methods.hh>

namespace Dumux
{
/*!
 * \ingroup SweModel
 *
 */
template <class TypeTag>
struct SweModelTraits
{
    using Indices = SweIndices;

    static constexpr int numEq() { return 3; }
    static constexpr int numPhases() { return 1; }
};


///////////////////////////////////////////////////////////////////////////
// properties for the shallow water model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the Swe inherits from shallow water
NEW_TYPE_TAG(Swe, INHERITS_FROM(ShallowWater));



//! Set friction law indices
NEW_PROP_TAG(Manning);
NEW_PROP_TAG(Chezy);
NEW_PROP_TAG(Nikuradse);
SET_INT_PROP(Swe, Manning, 1);
SET_INT_PROP(Swe, Chezy, 2);
SET_INT_PROP(Swe, Nikuradse, 3);

//! The model traits
SET_TYPE_PROP(Swe, ModelTraits, SweModelTraits<TypeTag>);

//! The local residual
SET_TYPE_PROP(Swe, LocalResidual, SweResidual<TypeTag>);

//! The volume variables
SET_TYPE_PROP(Swe, VolumeVariables, SweVolumeVariables<TypeTag>);

//! The flux variables
SET_TYPE_PROP(Swe, FluxVariables, SweFluxVariables<TypeTag>);

//! The flux variables cache
SET_TYPE_PROP(Swe, FluxVariablesCache, SweFluxVariablesCache<TypeTag>);

//! The specific vtk output fields
SET_TYPE_PROP(Swe, VtkOutputFields, SweVtkOutputFields);

//Set the advective Flux
SET_TYPE_PROP(Swe, AdvectionType, SweAdvectiveFlux<TypeTag>);
}

} // end namespace

#endif // DUMUX_SWE_MODEL_HH
