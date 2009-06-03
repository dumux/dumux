/*****************************************************************************
 *   Copyright (C) 2008 by Klaus Mosthaf, Andreas Lauser, Bernd Flemisch     *
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
#ifndef DUMUX_2P2CNIPROPERTIES_HH
#define DUMUX_2P2CNIPROPERTIES_HH

#include <dumux/boxmodels/2p2c/2p2cproperties.hh>

#include "2p2cnivertexdata.hh"
#include "2p2cnielementdata.hh"
#include "2p2cnifluxdata.hh"

namespace Dune
{
////////////////////////////////
// forward declarations
////////////////////////////////
template<class TypeTag>
class TwoPTwoCNIBoxModel;

template<class TypeTag>
class TwoPTwoCNIBoxJacobian;

template <class TypeTag>
class TwoPTwoCNIVertexData;

template <class TypeTag>
class TwoPTwoCNIElementData;

template <class TypeTag>
class TwoPTwoCNIFluxData;

/*!
 * \brief Enumerations for the non-isothermal 2-phase 2-component model
 */
template <int PVOffset = 0>
class TwoPTwoCNIIndices : public TwoPTwoCIndices<PVOffset>
{
public:
    static const int temperatureIdx = PVOffset + 2; //! The index for temperature in solution vectors.
};

////////////////////////////////
// properties
////////////////////////////////

namespace Properties
{
SET_INT_PROP(BoxTwoPTwoCNI, NumEq,         3); //!< set the number of equations to 3

//! Use the 2p2cni local jacobian operator for the 2p2cni model
SET_TYPE_PROP(BoxTwoPTwoCNI, 
              LocalJacobian,
              TwoPTwoCNIBoxJacobian<TypeTag>);

//! the Model property
SET_TYPE_PROP(BoxTwoPTwoCNI, Model, TwoPTwoCNIBoxModel<TypeTag>);

//! the VertexData property
SET_TYPE_PROP(BoxTwoPTwoCNI, VertexData, TwoPTwoCNIVertexData<TypeTag>);

//! the ElementData property
SET_TYPE_PROP(BoxTwoPTwoCNI, ElementData, TwoPTwoCNIElementData<TypeTag>);

//! the FluxData property
SET_TYPE_PROP(BoxTwoPTwoCNI, FluxData, TwoPTwoCNIFluxData<TypeTag>);

//! The indices required by the non-isothermal 2p2c model
SET_TYPE_PROP(BoxTwoPTwoCNI, TwoPTwoCIndices,   TwoPTwoCNIIndices<0>);
SET_TYPE_PROP(BoxTwoPTwoCNI, TwoPTwoCNIIndices, TwoPTwoCNIIndices<0>);

}

}
#endif
