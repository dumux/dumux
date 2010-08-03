// $Id: 2p2cniproperties.hh 3784 2010-06-24 13:43:57Z bernd $
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
/*!
 * \file
 *
 * \brief Defines the properties required for the non-isothermal two-phase,
 * two-component BOX model.
 */
#ifndef DUMUX_2P2CNIPROPERTIES_HH
#define DUMUX_2P2CNIPROPERTIES_HH

#include <dumux/boxmodels/2p2c/2p2cproperties.hh>

#include "2p2cnisecondaryvars.hh"

#include "2p2cnifluxvars.hh"

namespace Dumux
{
/*!
 * \addtogroup TwoPTwoCNIModel
 */
// \{
////////////////////////////////
// forward declarations
////////////////////////////////
template<class TypeTag>
class TwoPTwoCNIModel;

template<class TypeTag>
class TwoPTwoCNILocalResidual;

template <class TypeTag>
class TwoPTwoCNISecondaryVars;

template <class TypeTag>
class TwoPTwoCNIFluxVars;

/*!
 * \brief Enumerations for the non-isothermal 2-phase 2-component model
 */
template <class TypeTag, int formulation, int PVOffset>
class TwoPTwoCNIIndices : public TwoPTwoCIndices<TypeTag, formulation, PVOffset>
{
public:
    static const int temperatureIdx = PVOffset + 2; //! The index for temperature in primary variable vectors.
    static const int energyEqIdx = PVOffset + 2; //! The index for energy in equation vectors.
};

////////////////////////////////
// properties
////////////////////////////////

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the non-isothermal two-phase, two-component problems
NEW_TYPE_TAG(BoxTwoPTwoCNI, INHERITS_FROM(BoxTwoPTwoC));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(TwoPTwoCNIIndices); //!< Enumerations for the 2p2cni models

//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

SET_INT_PROP(BoxTwoPTwoCNI, NumEq,         3); //!< set the number of equations to 3

//! Use the 2p2cni local jacobian operator for the 2p2cni model
SET_TYPE_PROP(BoxTwoPTwoCNI,
              LocalResidual,
              TwoPTwoCNILocalResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(BoxTwoPTwoCNI, Model, TwoPTwoCNIModel<TypeTag>);

//! the SecondaryVars property
SET_TYPE_PROP(BoxTwoPTwoCNI, SecondaryVars, TwoPTwoCNISecondaryVars<TypeTag>);




//! the FluxVars property
SET_TYPE_PROP(BoxTwoPTwoCNI, FluxVars, TwoPTwoCNIFluxVars<TypeTag>);

//! The indices required by the non-isothermal 2p2c model
SET_PROP(BoxTwoPTwoCNI, TwoPTwoCIndices)
{ typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPTwoCNIIndices)) type; };

SET_PROP(BoxTwoPTwoCNI, TwoPTwoCNIIndices)
{ private:
    enum { formulation = GET_PROP_VALUE(TypeTag, PTAG(Formulation)) };
public:
    typedef TwoPTwoCNIIndices<TypeTag, formulation, 0> type;
};

}

}
#endif
