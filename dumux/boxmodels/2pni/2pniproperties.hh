// $Id: 2pniproperties.hh 3784 2010-06-24 13:43:57Z bernd $
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
/*!
 * \file
 *
 * \brief Defines the properties required for the non-isotherm two-phase BOX model.
 */

#ifndef DUMUX_2PNI_PROPERTIES_HH
#define DUMUX_2PNI_PROPERTIES_HH

#include <dumux/boxmodels/2p/2pproperties.hh>
#include "2pnivolumevariables.hh"

#include "2pnifluxvariables.hh"

namespace Dumux
{
/*!
 * \addtogroup TwoPNIBoxModel
 */
// \{

////////////////////////////////
// forward declarations
////////////////////////////////
template<class TypeTag>
class TwoPNIModel;

template<class TypeTag>
class TwoPNILocalResidual;

template <class TypeTag>
class TwoPNIVolumeVariables;

template <class TypeTag>
class TwoPNIFluxVariables;

/*!
 * \brief Enumerations for the non-isothermal two-phase model
 */
template <int PVOffset = 0>
class TwoPNIIndices : public TwoPIndices<PVOffset>
{
public:
    static const int temperatureIdx = PVOffset + 2; //! The primary variable index for temperature
    static const int energyEqIdx = PVOffset + 2; //! The equation index of the energy equation
};

////////////////////////////////
// properties
////////////////////////////////

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the non-isothermal two-phase problems
NEW_TYPE_TAG(BoxTwoPNI, INHERITS_FROM(BoxTwoP));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(TwoPNIIndices); //!< Enumerations for the non-isothermal 2p models

//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

SET_INT_PROP(BoxTwoPNI, NumEq, 3); //!< set the number of equations to 3

//! Use the 2pni local jacobian operator for the 2pni model
SET_TYPE_PROP(BoxTwoPNI,
              LocalResidual,
              TwoPNILocalResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(BoxTwoPNI, Model, TwoPNIModel<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(BoxTwoPNI, VolumeVariables, TwoPNIVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(BoxTwoPNI, FluxVariables, TwoPNIFluxVariables<TypeTag>);

//! The indices required by the non-isothermal two-phase model
SET_TYPE_PROP(BoxTwoPNI, TwoPIndices, TwoPNIIndices<0>);
SET_TYPE_PROP(BoxTwoPNI, TwoPNIIndices, TwoPNIIndices<0>);

}

}

#endif
