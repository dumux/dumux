// $Id:$
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
 * \brief Defines the properties required for the single-phase,
 *        two-compenent BOX model.
 */

#ifndef DUMUX_1P2C_PROPERTIES_HH
#define DUMUX_1P2C_PROPERTIES_HH

#include <dumux/boxmodels/tags.hh>

namespace Dune
{
////////////////////////////////
// forward declarations
////////////////////////////////
template<class TypeTag>
class OnePTwoCBoxModel;

template<class TypeTag>
class OnePTwoCBoxJacobian;

template <class TypeTag>
class OnePTwoCVertexData;

template <class TypeTag>
class OnePTwoCElementData;

template <class TypeTag>
class OnePTwoCFluxData;

/*!
 * \brief The indices for the isothermal single-phase, two-component model.
 */
struct OnePTwoCIndices
{
    // Primary variable indices
    static const int konti = 0;       //!< pressure in the solution vector
    static const int transport = 1;   //!< mole fraction in the solution vector
};

////////////////////////////////
// properties
////////////////////////////////
namespace Properties
{
SET_INT_PROP(BoxOnePTwoC, NumEq,         2); //!< set the number of equations to 2
SET_INT_PROP(BoxOnePTwoC, NumPhases,     1); //!< The number of phases in the 1p2c model is 1
SET_INT_PROP(BoxOnePTwoC, NumComponents, 2); //!< The number of components in the 1p2c model is 2

//! Set the default formulation to pWsN
SET_INT_PROP(BoxOnePTwoC, 
             Formulation,
             GET_PROP_TYPE(TypeTag,
                           PTAG(OnePTwoCIndices))::pWsN);

//! Use the 1p2c local jacobian operator for the 1p2c model
SET_TYPE_PROP(BoxOnePTwoC, 
              LocalJacobian,
              OnePTwoCBoxJacobian<TypeTag>);

//! the Model property
SET_TYPE_PROP(BoxOnePTwoC, Model, OnePTwoCBoxModel<TypeTag>);

//! the VertexData property
SET_TYPE_PROP(BoxOnePTwoC, VertexData, OnePTwoCVertexData<TypeTag>);

//! the ElementData property
SET_TYPE_PROP(BoxOnePTwoC, ElementData, OnePTwoCElementData<TypeTag>);

//! the FluxData property
SET_TYPE_PROP(BoxOnePTwoC, FluxData, OnePTwoCFluxData<TypeTag>);

//! the default upwind factor. Default 1.0, i.e. fully upwind...
SET_SCALAR_PROP(BoxOnePTwoC, UpwindAlpha, 1.0);

//! The indices required by the isothermal 1p2c model
SET_TYPE_PROP(BoxOnePTwoC, OnePTwoCIndices, Dune::OnePTwoCIndices);
}

}

#endif

