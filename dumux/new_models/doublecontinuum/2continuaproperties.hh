/*****************************************************************************
 *   Copyright (C) 2009 by Karin Erbertseder                                 *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Defines the properties required for the double continuum
 *        two-compenent BOX model.
 */

#ifndef DUMUX_TWO_CONTINUA_PROPERTIES_HH
#define DUMUX_TWO_CONTINUA_PROPERTIES_HH

#include <dumux/new_models/tags.hh>

namespace Dune
{
////////////////////////////////
// forward declarations
////////////////////////////////
template<class TypeTag>
class TwoContinuaBoxModel;

template<class TypeTag>
class TwoContinuaBoxJacobian;

template <class TypeTag>
class TwoContinuaVertexData;

template <class TypeTag>
class TwoContinuaElementData;

template <class TypeTag>
class TwoContinuaFluxData;

/*!
 * \brief The indices for the isothermal single-phase, two-component model.
 */
struct TwoContinuaIndices
{
    // Primary variable indices
    static const int pressure_bloodCont 	 = 0;
    static const int molefraction_bloodCont  = 1;
    static const int pressure_tissueCont 	 = 2;
    static const int molefraction_tissueCont = 3;

    // equation indices inside a continuum
    static const int konti = 0;
    static const int transport = 1;
    
    // indices of the continua
    static const int bloodCont = 0;	 //!< index of the blood continuum 
    static const int tissueCont = 1; //!< index of the tissue continuum
};

////////////////////////////////
// properties
////////////////////////////////
namespace Properties
{
SET_INT_PROP(BoxTwoContinua, NumEq,         4); //!< set the number of equations to 2
SET_INT_PROP(BoxTwoContinua, NumContinua,   2); //!< The number of phases in the double continuum model is 1
SET_INT_PROP(BoxTwoContinua, NumPhases,     1); //!< The number of phases in the double continuum model is 1
SET_INT_PROP(BoxTwoContinua, NumComponents, 2); //!< The number of components in the double continuum model is 2

//! Set the default formulation to pWsN
SET_INT_PROP(BoxTwoContinua, 
             Formulation,
             GET_PROP_TYPE(TypeTag,
                           PTAG(TwoContinuaIndices))::pWsN);

//! Use the double continuum local jacobian operator for the double continuum model
SET_TYPE_PROP(BoxTwoContinua, 
              LocalJacobian,
              TwoContinuaBoxJacobian<TypeTag>);

//! the Model property
SET_TYPE_PROP(BoxTwoContinua, Model, TwoContinuaBoxModel<TypeTag>);

//! the VertexData property
SET_TYPE_PROP(BoxTwoContinua, VertexData, TwoContinuaVertexData<TypeTag>);

//! the ElementData property
SET_TYPE_PROP(BoxTwoContinua, ElementData, TwoContinuaElementData<TypeTag>);

//! the FluxData property
SET_TYPE_PROP(BoxTwoContinua, FluxData, TwoContinuaFluxData<TypeTag>);

//! the default upwind factor. Default 1.0, i.e. fully upwind...
SET_SCALAR_PROP(BoxTwoContinua, UpwindAlpha, 1.0);

//! The indices required by the isothermal double continuum model
SET_TYPE_PROP(BoxTwoContinua, TwoContinuaIndices, Dune::TwoContinuaIndices);
}

}

#endif

