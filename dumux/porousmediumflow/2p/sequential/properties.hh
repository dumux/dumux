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
 * \ingroup SequentialTwoPModel
 * \brief Defines the properties required for (immiscible) two-phase sequential models.
 */
#ifndef DUMUX_2PPROPERTIES_SEQUENTIAL_HH
#define DUMUX_2PPROPERTIES_SEQUENTIAL_HH

//Dumux-includes
#include <dumux/porousmediumflow/sequential/properties.hh>
#include "indices.hh"
#include <dumux/material/spatialparams/sequentialfv.hh>

namespace Dumux
{
////////////////////////////////
// forward declarations
////////////////////////////////
template <class TypeTag, bool enableCompressibility>
class CellData2P;

////////////////////////////////
// properties
////////////////////////////////
namespace Properties
{
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The TypeTag for sequential two-phase problems
NEW_TYPE_TAG(SequentialTwoP, INHERITS_FROM(SequentialModel));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
NEW_PROP_TAG( SpatialParams ); //!< The type of the spatial parameters object
NEW_PROP_TAG(MaterialLaw);   //!< The material law which ought to be used (extracted from the spatial parameters)
NEW_PROP_TAG( Formulation); //!< The formulation of the model
NEW_PROP_TAG( PressureFormulation); //!< The formulation of the pressure model
NEW_PROP_TAG( SaturationFormulation); //!< The formulation of the saturation model
NEW_PROP_TAG( VelocityFormulation); //!< The type of velocity reconstructed for the transport model
NEW_PROP_TAG( EnableCompressibility);//!< Returns whether compressibility is allowed
NEW_PROP_TAG( FluidSystem ); //!< Defines the fluid system
NEW_PROP_TAG( FluidState );//!< Defines the fluid state
}
}

#include <dumux/porousmediumflow/sequential/variableclass.hh>
#include <dumux/porousmediumflow/2p/sequential/celldata.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>
#include <dumux/material/fluidstates/isothermalimmiscible.hh>

namespace Dumux
{
namespace Properties
{
//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////
//! Set number of equations to 2 for isothermal two-phase models
SET_INT_PROP(SequentialTwoP, NumEq, 2);

//! Set number of phases to 2 for two-phase models
SET_INT_PROP(SequentialTwoP, NumPhases, 2);//!< The number of phases in the 2p model is 2

//! Set number of components to 1 for immiscible two-phase models
SET_INT_PROP(SequentialTwoP, NumComponents, 1); //!< Each phase consists of 1 pure component

//! Set \f$p_w\f$-\f$S_w\f$ formulation as default two-phase formulation
SET_INT_PROP(SequentialTwoP, Formulation, SequentialTwoPCommonIndices::pwsw);

//! Chose the set of indices depending on the chosen formulation
SET_PROP(SequentialTwoP, Indices)
{
    using type = SequentialTwoPIndices<GET_PROP_VALUE(TypeTag, Formulation), 0>;
};

//! Set the default pressure formulation according to the chosen two-phase formulation
SET_INT_PROP(SequentialTwoP,
    PressureFormulation,
    GET_PROP_TYPE(TypeTag, Indices)::pressureType);

//! Set the default saturation formulation according to the chosen two-phase formulation
SET_INT_PROP(SequentialTwoP,
    SaturationFormulation,
    GET_PROP_TYPE(TypeTag, Indices)::saturationType);

//! Set the default velocity formulation according to the chosen two-phase formulation
SET_INT_PROP(SequentialTwoP,
    VelocityFormulation,
    GET_PROP_TYPE(TypeTag, Indices)::velocityDefault);

//! Disable compressibility by default
SET_BOOL_PROP(SequentialTwoP, EnableCompressibility, false);

//! Set general sequential VariableClass as default
SET_TYPE_PROP(SequentialTwoP, Variables, VariableClass<TypeTag>);

//! Set standart CellData of immiscible two-phase models as default
SET_TYPE_PROP(SequentialTwoP, CellData, CellData2P<TypeTag, GET_PROP_VALUE(TypeTag, EnableCompressibility)>);

//! Set default fluid state
SET_PROP(SequentialTwoP, FluidState)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
public:
    using type = IsothermalImmiscibleFluidState<Scalar, FluidSystem>;
};

//! The spatial parameters to be employed. Use SequentialFVSpatialParams by default.
SET_TYPE_PROP(SequentialTwoP, SpatialParams, SequentialFVSpatialParams<TypeTag>);
// \}
}

}

#endif
