// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Markus Wolff                                 *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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
 * \ingroup IMPES
 * \ingroup IMPETProperties
 */
/*!
 * \file
 *
 * \brief Defines the properties required for (immiscible) two-phase sequential models.
 */

#ifndef DUMUX_2PPROPERTIES_DECOUPLED_HH
#define DUMUX_2PPROPERTIES_DECOUPLED_HH

//Dumux-includes
#include <dumux/decoupled/common/decoupledproperties.hh>
#include "2pindices.hh"

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

//! The TypeTag for decoupled two-phase problems
NEW_TYPE_TAG(DecoupledTwoP, INHERITS_FROM(DecoupledModel))
;

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG( SpatialParameters )
; //!< The type of the spatial parameters object
NEW_PROP_TAG(MaterialLaw);   //!< The material law which ought to be used (extracted from the spatial parameters)
NEW_PROP_TAG(MaterialLawParams); //!< The material law parameters (extracted from the spatial parameters)
NEW_PROP_TAG( EnableGravity)
; //!< Returns whether gravity is considered in the problem
NEW_PROP_TAG( Formulation); //!< The formulation of the model
NEW_PROP_TAG( PressureFormulation)
; //!< The formulation of the pressure model
NEW_PROP_TAG( SaturationFormulation)
; //!< The formulation of the saturation model
NEW_PROP_TAG( VelocityFormulation)
; //!< The type of velocity reconstructed for the transport model
NEW_PROP_TAG( EnableCompressibility)
;// !< Returns whether compressibility is allowed
NEW_PROP_TAG( WettingPhase)
; //!< The wetting phase of a two-phase model
NEW_PROP_TAG( NonwettingPhase)
; //!< The non-wetting phase of a two-phase model
NEW_PROP_TAG( FluidSystem )//!< Defines the fluid system
;
NEW_PROP_TAG( FluidState )//!< Defines the fluid state
;

//! \cond \private
// keep only for compatibility with box models
NEW_PROP_TAG( TwoPIndices );
//! \endcond

NEW_PROP_TAG( ErrorTermFactor ); //! Scaling factor for the error term (term to damp unphysical saturation overshoots via pressure correction)
NEW_PROP_TAG( ErrorTermLowerBound );//! Lower threshold used for the error term evaluation (term to damp unphysical saturation overshoots via pressure correction)
NEW_PROP_TAG( ErrorTermUpperBound );//!Upper threshold used for the error term evaluation (term to damp unphysical saturation overshoots via pressure correction)
}
}

#include <dumux/decoupled/common/variableclass.hh>
#include <dumux/decoupled/2p/cellData2p.hh>
#include <dumux/material/fluidsystems/2pimmisciblefluidsystem.hh>
#include <dumux/material/fluidstates/isothermalimmisciblefluidstate.hh>

namespace Dumux
{
namespace Properties
{
//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////
//! Set number of equations to 2 for isothermal two-phase models
SET_INT_PROP(DecoupledTwoP, NumEq, 2);

//! Set number of phases to 2 for two-phase models
SET_INT_PROP(DecoupledTwoP, NumPhases, 2);//!< The number of phases in the 2p model is 2

//! Set number of components to 1 for immiscible two-phase models
SET_INT_PROP(DecoupledTwoP, NumComponents, 1); //!< Each phase consists of 1 pure component

//! Set \f$p_w\f$-\f$S_n\f$ formulation as default two-phase formulation
SET_INT_PROP(DecoupledTwoP,
    Formulation,
    DecoupledTwoPCommonIndices::pwSw);

//! Chose the set of indices depending on the chosen formulation
SET_PROP(DecoupledTwoP, Indices)
{
typedef DecoupledTwoPIndices<GET_PROP_VALUE(TypeTag, Formulation), 0> type;
};

//! \cond \private
// keep only for compatibility with box models
SET_TYPE_PROP(DecoupledTwoP, TwoPIndices, typename GET_PROP_TYPE(TypeTag, Indices));
//! \endcond

//! Set the default pressure formulation according to the chosen two-phase formulation
SET_INT_PROP(DecoupledTwoP,
    PressureFormulation,
    GET_PROP_TYPE(TypeTag, Indices)::pressureType);

//! Set the default saturation formulation according to the chosen two-phase formulation
SET_INT_PROP(DecoupledTwoP,
    SaturationFormulation,
    GET_PROP_TYPE(TypeTag, Indices)::saturationType);

//! Set the default velocity formulation according to the chosen two-phase formulation
SET_INT_PROP(DecoupledTwoP,
    VelocityFormulation,
    GET_PROP_TYPE(TypeTag, Indices)::velocityDefault);

//! Disable compressibility by default
SET_BOOL_PROP(DecoupledTwoP, EnableCompressibility, false);

//! Set general decoupled VariableClass as default
SET_TYPE_PROP(DecoupledTwoP, Variables, VariableClass<TypeTag>);

//! Set standart CellData of immiscible two-phase models as default
SET_TYPE_PROP(DecoupledTwoP, CellData, CellData2P<TypeTag, GET_PROP_VALUE(TypeTag, EnableCompressibility)>);

//! Set default fluid system
SET_TYPE_PROP(DecoupledTwoP, FluidSystem, TwoPImmiscibleFluidSystem<TypeTag>);

//! Set default fluid state
SET_PROP(DecoupledTwoP, FluidState)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
public:
    typedef IsothermalImmiscibleFluidState<Scalar, FluidSystem> type;
};

/*!
 * \brief Set the property for the material parameters by extracting
 *        it from the material law.
 */
SET_PROP(DecoupledTwoP, MaterialLawParams)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;

public:
    typedef typename MaterialLaw::Params type;
};

//! Default error term factor
SET_SCALAR_PROP(DecoupledTwoP, ErrorTermFactor, 0.5);
//! Default lower threshold for evaluation of an error term
SET_SCALAR_PROP(DecoupledTwoP, ErrorTermLowerBound, 0.1);
//! Default upper threshold for evaluation of an error term
SET_SCALAR_PROP(DecoupledTwoP, ErrorTermUpperBound, 0.9);

// \}
}

}

#endif
