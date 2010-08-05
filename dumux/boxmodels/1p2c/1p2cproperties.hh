// $Id: 1p2cproperties.hh 3838 2010-07-15 08:31:53Z bernd $
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

#include<dumux/boxmodels/common/boxproperties.hh>

namespace Dumux
{
////////////////////////////////
// forward declarations
////////////////////////////////
template<class TypeTag>
class OnePTwoCBoxModel;

template<class TypeTag>
class OnePTwoCLocalResidual;

template <class TypeTag>
class OnePTwoCVolumeVariables;

template <class TypeTag>
class OnePTwoCFluxVariables;

/*!
 * \brief The indices for the isothermal single-phase, two-component model.
 */
struct OnePTwoCIndices
{
    // Primary variable indices
    static const int konti = 0;       //!< pressure in the solution vector
    static const int transport = 1;   //!< mole fraction in the solution vector
};

namespace Properties
{

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the isothermal single-phase, two-component problems
NEW_TYPE_TAG(BoxOnePTwoC, INHERITS_FROM(BoxModel));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(NumPhases);   //!< Number of fluid phases in the system
NEW_PROP_TAG(NumComponents);   //!< Number of fluid components in the system
NEW_PROP_TAG(OnePTwoCIndices); //!< Enumerations for the 1p2c models
NEW_PROP_TAG(SpatialParameters); //!< The type of the spatial parameters
NEW_PROP_TAG(FluidSystem); //!< Type of the multi-component relations
NEW_PROP_TAG(UpwindAlpha);   //!< The default value of the upwind parameter
NEW_PROP_TAG(EnableGravity); //!< Returns whether gravity is considered in the problem

//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

SET_INT_PROP(BoxOnePTwoC, NumEq, 2); //!< set the number of equations to 2
SET_INT_PROP(BoxOnePTwoC, NumPhases, 1); //!< The number of phases in the 1p2c model is 1
SET_INT_PROP(BoxOnePTwoC, NumComponents, 2); //!< The number of components in the 1p2c model is 2

//! Use the 1p2c local jacobian operator for the 1p2c model
SET_TYPE_PROP(BoxOnePTwoC,
              LocalResidual,
              OnePTwoCLocalResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(BoxOnePTwoC, Model, OnePTwoCBoxModel<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(BoxOnePTwoC, VolumeVariables, OnePTwoCVolumeVariables<TypeTag>);




//! the FluxVariables property
SET_TYPE_PROP(BoxOnePTwoC, FluxVariables, OnePTwoCFluxVariables<TypeTag>);

//! the default upwind factor. Default 1.0, i.e. fully upwind...
SET_SCALAR_PROP(BoxOnePTwoC, UpwindAlpha, 1.0);

//! The indices required by the isothermal 1p2c model
SET_TYPE_PROP(BoxOnePTwoC, OnePTwoCIndices, Dumux::OnePTwoCIndices);
}

}

#endif

