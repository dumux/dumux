/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
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
 * \brief Contains the properties for the Richards BOX model.
 */
#ifndef DUMUX_RICHARDS_PROPERTIES_DATA_HH
#define DUMUX_RICHARDS_PROPERTIES_DATA_HH

#include <dumux/boxmodels/tags.hh>

namespace Dune
{
/*!
 * \addtogroup RichardsBoxModel
 */
// \{
////////////////////////////////
// forward declarations
////////////////////////////////
template<class TypeTag>
class RichardsBoxModel;

template<class TypeTag>
class RichardsBoxJacobian;

template <class TypeTag>
class RichardsVertexData;

template <class TypeTag>
class RichardsElementData;

template <class TypeTag>
class RichardsFluxData;

/*!
 * \brief Indices for the single phase model.
 */
struct RichardsIndices
{
    static const int pWIdx = 0;
};

///////////////////////////////////////////////////////////////////////////
// properties for the isothermal single phase model
///////////////////////////////////////////////////////////////////////////
namespace Properties {
SET_INT_PROP(BoxRichards, NumEq, 1);
SET_INT_PROP(BoxRichards, NumPhases, 2);

//! Use the 2p local jacobian operator for the 2p model
SET_TYPE_PROP(BoxRichards,
              LocalJacobian,
              RichardsBoxJacobian<TypeTag>);

//! the Model property
SET_TYPE_PROP(BoxRichards, Model, RichardsBoxModel<TypeTag>);

//! the VertexData property
SET_TYPE_PROP(BoxRichards, VertexData, RichardsVertexData<TypeTag>);

//! the ElementData property
SET_TYPE_PROP(BoxRichards, ElementData, RichardsElementData<TypeTag>);

//! the FluxData property
SET_TYPE_PROP(BoxRichards, FluxData, RichardsFluxData<TypeTag>);

//! the default upwind factor. Default 1.0, i.e. fully upwind...
SET_SCALAR_PROP(BoxRichards, UpwindAlpha, 1.0);

//! the upwind factor for the mobility. uses the value of UpwindAlpha
//! if the property is not overwritten elsewhere
SET_SCALAR_PROP(BoxRichards,
                MobilityUpwindAlpha,
                GET_PROP_VALUE(TypeTag, PTAG(UpwindAlpha)));

//! The indices required by the isothermal single-phase model
SET_TYPE_PROP(BoxRichards, RichardsIndices, Dune::RichardsIndices);

// \}
};

} // end namepace

#endif
