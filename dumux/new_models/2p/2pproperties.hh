//$Id:$
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
 * \brief Defines the properties required for the twophase BOX model.
 */

#ifndef DUMUX_2PPROPERTIES_HH
#define DUMUX_2PPROPERTIES_HH

#include <dumux/new_models/tags.hh>

namespace Dune
{
////////////////////////////////
// forward declarations
////////////////////////////////
template<class TypeTag>
class TwoPBoxModel;

template<class TypeTag>
class TwoPBoxJacobian;

template <class TypeTag>
class TwoPVertexData;

template <class TypeTag>
class TwoPElementData;

template <class TypeTag>
class TwoPFluxData;

/*!
 * \brief The indices for the isothermal two-phase model.
 *
 * \tparam PVOffset    The first index in a primary variable vector.
 */
template <int PVOffset = 0>
struct TwoPIndices
{
    // Primary variable indices
    static const int pressureIdx   = PVOffset + 0; //!< Index for wetting/non-wetting phase pressure (depending on formulation) in a solution vector
    static const int saturationIdx = PVOffset + 1; //!< Index of the saturation of the non-wetting/wetting phase
   
    // Formulations
    static const int pWsN = 0; //!< Pw and Sn as primary variables
    static const int pNsW = 1; //!< Pn and Sw as primary variables

    // Phase indices
    static const int wPhase = 0; //!< Index of the wetting phase in a phase vector
    static const int nPhase = 1; //!< Index of the non-wetting phase in a phase vector

    /*!
     * \brief Convert a index in a phase vector to an index in a
     *        vector of primary variables.
     */
    static int phase2Mass(int phaseIdx)
    {
        return PVOffset + phaseIdx;
    };
};

////////////////////////////////
// properties
////////////////////////////////
namespace Properties
{
SET_INT_PROP(BoxTwoP, NumEq,         2); //!< set the number of equations to 2
SET_INT_PROP(BoxTwoP, NumPhases,     2); //!< The number of phases in the 2p model is 2

//! Set the default formulation to pWsN
SET_INT_PROP(BoxTwoP, 
             Formulation,
             GET_PROP_TYPE(TypeTag,
                           PTAG(TwoPIndices))::pWsN);

//! Use the 2p local jacobian operator for the 2p model
SET_TYPE_PROP(BoxTwoP, 
              LocalJacobian,
              TwoPBoxJacobian<TypeTag>);

//! the Model property
SET_TYPE_PROP(BoxTwoP, Model, TwoPBoxModel<TypeTag>);

//! the VertexData property
SET_TYPE_PROP(BoxTwoP, VertexData, TwoPVertexData<TypeTag>);

//! the ElementData property
SET_TYPE_PROP(BoxTwoP, ElementData, TwoPElementData<TypeTag>);

//! the FluxData property
SET_TYPE_PROP(BoxTwoP, FluxData, TwoPFluxData<TypeTag>);

//! the default upwind factor. Default 1.0, i.e. fully upwind...
SET_SCALAR_PROP(BoxTwoP, UpwindAlpha, 1.0);

//! the upwind factor for the mobility. uses the value of UpwindAlpha
//! if the property is not overwritten elsewhere
SET_SCALAR_PROP(BoxTwoP, 
                MobilityUpwindAlpha, 
                GET_PROP_VALUE(TypeTag, PTAG(UpwindAlpha)));

//! The indices required by the isothermal 2p model
SET_TYPE_PROP(BoxTwoP, TwoPIndices, TwoPIndices<0>);
}

}

#endif

