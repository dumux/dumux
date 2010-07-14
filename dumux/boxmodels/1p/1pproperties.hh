// $Id: 1pproperties.hh 3784 2010-06-24 13:43:57Z bernd $
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
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
 * \brief Defines the properties required for the onephase BOX model.
 */
#ifndef DUMUX_1P_PROPERTIES_DATA_HH
#define DUMUX_1P_PROPERTIES_DATA_HH

#include <dumux/boxmodels/boxscheme/boxproperties.hh>

namespace Dumux
{
/*!
 * \addtogroup OnePBoxModel
 */
// \{
////////////////////////////////
// forward declarations
////////////////////////////////
template<class TypeTag>
class OnePBoxModel;

template<class TypeTag>
class OnePBoxJacobian;

template <class TypeTag>
class OnePVertexData;

template <class TypeTag>
class OnePElementData;

template <class TypeTag>
class OnePFluxData;

/*!
 * \brief Indices for the single phase model.
 */
struct OnePIndices
{
    static const int pressureIdx = 0;
};

///////////////////////////////////////////////////////////////////////////
// properties for the isothermal single phase model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the isothermal single phase problems
NEW_TYPE_TAG(BoxOneP, INHERITS_FROM(BoxScheme));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(NumPhases);   //!< Number of fluid phases in the system
NEW_PROP_TAG(OnePIndices); //!< Enumerations for the 1p models
NEW_PROP_TAG(SpatialParameters); //!< The type of the soil properties object
NEW_PROP_TAG(Fluid); //!< The fluid for the single-phase problems
NEW_PROP_TAG(EnableGravity); //!< Returns whether gravity is considered in the problem

//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

SET_INT_PROP(BoxOneP, NumEq, 1);
SET_INT_PROP(BoxOneP, NumPhases, 1);

//! Use the 2p local jacobian operator for the 2p model
SET_TYPE_PROP(BoxOneP,
              LocalJacobian,
              OnePBoxJacobian<TypeTag>);

//! the Model property
SET_TYPE_PROP(BoxOneP, Model, OnePBoxModel<TypeTag>);

//! the VertexData property
SET_TYPE_PROP(BoxOneP, VertexData, OnePVertexData<TypeTag>);

//! the ElementData property
SET_TYPE_PROP(BoxOneP, ElementData, OnePElementData<TypeTag>);

//! the FluxData property
SET_TYPE_PROP(BoxOneP, FluxData, OnePFluxData<TypeTag>);

//! The indices required by the isothermal single-phase model
SET_TYPE_PROP(BoxOneP, OnePIndices, OnePIndices);

// \}
};

} // end namepace

#endif
