//$Id$
/*****************************************************************************
 *   Copyright (C) 2008 by Klaus Mosthaf                                     *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Copyright (C) 2008 Bernd Flemisch                                       *
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

#ifndef DUMUX_2PNI_PROPERTIES_HH
#define DUMUX_2PNI_PROPERTIES_HH

#include <dumux/new_models/2p/2pproperties.hh>

#include "2pnivertexdata.hh"
#include "2pnielementdata.hh"
#include "2pnifluxdata.hh"

namespace Dune
{
////////////////////////////////
// forward declarations
////////////////////////////////
template<class TypeTag>
class TwoPNIBoxModel;

template<class TypeTag>
class TwoPNIBoxJacobian;

template <class TypeTag>
class TwoPNIVertexData;

template <class TypeTag>
class TwoPNIElementData;

template <class TypeTag>
class TwoPNIFluxData;

/*!
 * \brief Enumerations for the non-isothermal 2-phase model
 */
template <int PVOffset = 0>
class TwoPNIIndices : public TwoPIndices<PVOffset>
{
public:
    static const int temperatureIdx = PVOffset + 2; //! The index for temperature in solution vectors.
};

////////////////////////////////
// properties
////////////////////////////////

namespace Properties
{
SET_INT_PROP(BoxTwoPNI, NumEq,         3); //!< set the number of equations to 3

//! Use the 2p2cni local jacobian operator for the 2p2cni model
SET_TYPE_PROP(BoxTwoPNI, 
              LocalJacobian,
              TwoPNIBoxJacobian<TypeTag>);

//! the Model property
SET_TYPE_PROP(BoxTwoPNI, Model, TwoPNIBoxModel<TypeTag>);

//! the VertexData property
SET_TYPE_PROP(BoxTwoPNI, VertexData, TwoPNIVertexData<TypeTag>);

//! the ElementData property
SET_TYPE_PROP(BoxTwoPNI, ElementData, TwoPNIElementData<TypeTag>);

//! the FluxData property
SET_TYPE_PROP(BoxTwoPNI, FluxData, TwoPNIFluxData<TypeTag>);

//! The indices required by the non-isothermal 2p2c model
SET_TYPE_PROP(BoxTwoPNI, TwoPIndices,   TwoPNIIndices<0>);
SET_TYPE_PROP(BoxTwoPNI, TwoPNIIndices, TwoPNIIndices<0>);

}

}

#endif
