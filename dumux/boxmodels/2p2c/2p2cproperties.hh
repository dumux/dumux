// $Id$
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

#ifndef DUMUX_2P2CTRAITS_HH
#define DUMUX_2P2CTRAITS_HH

#include "2p2cnewtoncontroller.hh"

#include "2p2cvertexdata.hh"
#include "2p2celementdata.hh"
#include "2p2cfluxdata.hh"

namespace Dune
{
////////////////////////////////
// forward declarations
////////////////////////////////
template<class TypeTag>
class TwoPTwoCBoxModel;

template<class TypeTag>
class TwoPTwoCBoxJacobian;

template <class TypeTag>
class TwoPTwoCVertexData;

template <class TypeTag>
class TwoPTwoCElementData;

template <class TypeTag>
class TwoPTwoCFluxData;

/*!
 * \brief The formulation independent indices for the 2p2c model which
 *        do not depend on an offset in the primary variable vector.
 */
struct TwoPTwoCCommonIndices
{
    // Phase state (-> 'pseudo' primary variable)
    static const int nPhaseOnly = 0; //!< Only the non-wetting phase is present
    static const int wPhaseOnly = 1; //!< Only the wetting phase is present
    static const int bothPhases = 2; //!< Both phases are present
    
    // Formulations
    static const int pWsN = 0; //!< Pw and Sn as primary variables
    static const int pNsW = 1; //!< Pn and Sw as primary variables

    // Phase indices
    static const int wPhase = 0; //!< Index of the wetting phase in a phase vector
    static const int nPhase = 1; //!< Index of the non-wetting phase in a phase vector

    // Component indices
    static const int wComp = 0; //!< Index of the wetting component in a component vector
    static const int nComp = 1; //!< Index of the non-wetting component in a compent vector
};

/*!
 * \brief The indices for the isothermal TwoPTwoC model.
 *
 * \tparam PVOffset    The first index in a primary variable vector.
 */
template <int formulation = TwoPTwoCCommonIndices::pWsN, int PVOffset = 0>
struct TwoPTwoCIndices
    : public TwoPTwoCCommonIndices
{
    // Primary variable indices
    static const int pressureIdx = PVOffset + 0; //!< Index for wetting/non-wetting phase pressure (depending on formulation) in a solution vector
    static const int switchIdx   = PVOffset + 1; //!< Index of the either the saturation or the mass fraction of the non-wetting/wetting phase

    static const int pW = pressureIdx; //!< Index for wetting/non-wetting phase pressure (depending on formulation) in a solution vector
    static const int sNorX = switchIdx; //!< Index of the either the saturation or the mass fraction of the non-wetting/wetting phase
  
    /*!
     * \brief Map a component index to a mass index.
     *
     * (The mass index is the index of a component in the result
     * vector of primary variables in the storage or flux terms.)
     */
    static int comp2Mass(int compIdx) { return PVOffset + compIdx; }
};

/*!
 * \brief The indices for the isothermal TwoPTwoC model.
 *
 * \tparam PVOffset    The first index in a primary variable vector.
 */
template <int PVOffset>
struct TwoPTwoCIndices<TwoPTwoCCommonIndices::pNsW, PVOffset> 
    : public TwoPTwoCCommonIndices
{
    // Primary variable indices
    static const int pressureIdx = PVOffset + 0; //!< Index for wetting/non-wetting phase pressure (depending on formulation) in a solution vector
    static const int switchIdx   = PVOffset + 1; //!< Index of the either the saturation or the mass fraction of the non-wetting/wetting phase

    static const int pN = pressureIdx; //!< Index for wetting/non-wetting phase pressure (depending on formulation) in a solution vector
    static const int sWorX = switchIdx; //!< Index of the either the saturation or the mass fraction of the non-wetting/wetting phase
  
    /*!
     * \brief Map a component index to a mass index.
     *
     * (The mass index is the index of a component in the result
     * vector of primary variables in the storage or flux terms.)
     */
    static int comp2Mass(int compIdx) { return PVOffset + 1 - compIdx; }
};

////////////////////////////////
// properties
////////////////////////////////

namespace Properties
{
SET_INT_PROP(BoxTwoPTwoC, NumEq,         2); //!< set the number of equations to 2
SET_INT_PROP(BoxTwoPTwoC, NumPhases,     2); //!< The number of phases in the 2p2c model is 2
SET_INT_PROP(BoxTwoPTwoC, NumComponents, 2); //!< The number of components in the 2p2c model is 2

//! Set the default formulation to pWsN
SET_INT_PROP(BoxTwoPTwoC, 
             Formulation,
             TwoPTwoCCommonIndices::pWsN);

//! Use the 2p2c local jacobian operator for the 2p2c model
SET_TYPE_PROP(BoxTwoPTwoC, 
              LocalJacobian,
              TwoPTwoCBoxJacobian<TypeTag>);

//! Use the 2p2c specific newton controller for the 2p2c model
SET_PROP(BoxTwoPTwoC, NewtonController)
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonMethod))  NewtonMethod;
  
public:
    typedef TwoPTwoCNewtonController<NewtonMethod, TypeTag> type;
};

//! the Model property
SET_TYPE_PROP(BoxTwoPTwoC, Model, TwoPTwoCBoxModel<TypeTag>);

//! the VertexData property
SET_TYPE_PROP(BoxTwoPTwoC, VertexData, TwoPTwoCVertexData<TypeTag>);

//! the ElementData property
SET_TYPE_PROP(BoxTwoPTwoC, ElementData, TwoPTwoCElementData<TypeTag>);

//! the FluxData property
SET_TYPE_PROP(BoxTwoPTwoC, FluxData, TwoPTwoCFluxData<TypeTag>);

//! the default upwind factor. Default 1.0, i.e. fully upwind...
SET_SCALAR_PROP(BoxTwoPTwoC, UpwindAlpha, 1.0);

//! the upwind factor for the mobility. uses the value of UpwindAlpha
//! if the property is not overwritten elsewhere
SET_SCALAR_PROP(BoxTwoPTwoC, 
                MobilityUpwindAlpha, 
                GET_PROP_VALUE(TypeTag, PTAG(UpwindAlpha)));

//! The indices required by the isothermal 2p2c model
SET_PROP(BoxTwoPTwoC, 
         TwoPTwoCIndices)
{
    typedef TwoPTwoCIndices<GET_PROP_VALUE(TypeTag, PTAG(Formulation)), 0> type;
};

}

}

#endif
