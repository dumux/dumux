// $Id: 2pproperties.hh 3357 2010-03-25 13:02:05Z lauser $
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

namespace Dumux
{
/*!
 * \addtogroup TwoPBoxModel
 */
// \{

////////////////////////////////
// forward declarations
////////////////////////////////
template<class TypeTag>
class TwoPModel;

template<class TypeTag>
class TwoPLocalResidual;

template<class TypeTag>
class TwoPNewtonController;

template <class TypeTag>
class TwoPProblem;

template <class TypeTag>
class TwoPVolumeVariables;

template <class TypeTag>
class TwoPFluxVariables;

template<class TypeTag>
class FluidSystem2P;

template<class TypeTag>
class TwoPFluidState;

/*!
 * \brief The common indices for the isothermal two-phase model.
 */
struct TwoPCommonIndices
{
    // Formulations
    static const int pwSn = 0; //!< Pw and Sn as primary variables
    static const int pnSw = 1; //!< Pn and Sw as primary variables

    // Phase indices
    static const int wPhaseIdx = 0; //!< Index of the wetting phase in a phase vector
    static const int nPhaseIdx = 1; //!< Index of the non-wetting phase in a phase vector
};

/*!
 * \brief The indices for the \f$p_w-S_n\f$ formulation of the
 *        isothermal two-phase model.
 *
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <int formulation = TwoPCommonIndices::pwSn, int PVOffset = 0>
struct TwoPIndices : public TwoPCommonIndices
{
    // Primary variable indices
    static const int pressureIdx = PVOffset + 0; //!< Index for wetting/non-wetting phase pressure (depending on formulation) in a solution vector
    static const int saturationIdx = PVOffset + 1; //!< Index of the saturation of the non-wetting/wetting phase

    // indices of the primary variables
    static const int pwIdx = PVOffset + 0; //!< Pressure index of the wetting phase
    static const int SnIdx = PVOffset + 1; //!< Saturation index of the wetting phase

    // indices of the equations
    static const int contiWEqIdx = PVOffset + 0; //!< Index of the continuity equation of the wetting phase
    static const int contiNEqIdx = PVOffset + 1; //!< Index of the continuity equation of the non-wetting phase
};

/*!
 * \brief The indices for the \f$p_w-S_n\f$ formulation of the
 *        isothermal two-phase model.
 *
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <int PVOffset>
struct TwoPIndices<TwoPCommonIndices::pnSw, PVOffset>
    : public TwoPCommonIndices
{
    // Primary variable indices
    static const int pressureIdx = PVOffset + 0; //!< Index for wetting/non-wetting phase pressure (depending on formulation) in a solution vector
    static const int saturationIdx = PVOffset + 1; //!< Index of the saturation of the non-wetting/wetting phase

    // indices of the primary variables
    static const int pnIdx = PVOffset + 0; //!< Pressure index of the wetting phase
    static const int SwIdx = PVOffset + 1; //!< Saturation index of the wetting phase

    // indices of the equations
    static const int contiNEqIdx = PVOffset + 0; //!< Index of the continuity equation of the non-wetting phase
    static const int contiWEqIdx = PVOffset + 1; //!< Index of the continuity equation of the wetting phase
};

// \}

////////////////////////////////
// properties
////////////////////////////////
namespace Properties
{

/*!
 * \addtogroup TwoPBoxModel
 */
// \{

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the two-phase problems
NEW_TYPE_TAG(BoxTwoP, INHERITS_FROM(BoxModel));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(NumPhases);   //!< Number of fluid phases in the system
NEW_PROP_TAG(EnableGravity); //!< Returns whether gravity is considered in the problem
NEW_PROP_TAG(MobilityUpwindAlpha); //!< The value of the upwind parameter for the mobility
NEW_PROP_TAG(Formulation);   //!< The formulation of the model
NEW_PROP_TAG(TwoPIndices); //!< Enumerations for the 2p models
NEW_PROP_TAG(SpatialParameters); //!< The type of the spatial parameters object
NEW_PROP_TAG(MaterialLaw);   //!< The material law which ought to be used (extracted from the spatial parameters)
NEW_PROP_TAG(MaterialLawParams); //!< The context material law (extracted from the spatial parameters)
NEW_PROP_TAG(WettingPhase); //!< The wetting phase for two-phase models
NEW_PROP_TAG(NonwettingPhase); //!< The non-wetting phase for two-phase models
NEW_PROP_TAG( FluidSystem ); //!<The fluid systems including the information about the phases
NEW_PROP_TAG( FluidState ); //!<The phases state

//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

SET_INT_PROP(BoxTwoP, NumEq, 2); //!< set the number of equations to 2
SET_INT_PROP(BoxTwoP, NumPhases, 2); //!< The number of phases in the 2p model is 2

//! Set the default formulation to pWsN
SET_INT_PROP(BoxTwoP,
             Formulation,
             TwoPCommonIndices::pwSn);

//! Use the 2p local jacobian operator for the 2p model
SET_TYPE_PROP(BoxTwoP,
              LocalResidual,
              TwoPLocalResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(BoxTwoP, Model, TwoPModel<TypeTag>);

//! the default newton controller for two-phase problems
SET_TYPE_PROP(BoxTwoP, NewtonController, TwoPNewtonController<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(BoxTwoP, VolumeVariables, TwoPVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(BoxTwoP, FluxVariables, TwoPFluxVariables<TypeTag>);

//! the upwind factor for the mobility.
SET_SCALAR_PROP(BoxTwoP, MobilityUpwindAlpha, 1.0);

//! The indices required by the isothermal 2p model
SET_PROP(BoxTwoP, TwoPIndices)
{
    typedef TwoPIndices<GET_PROP_VALUE(TypeTag, PTAG(Formulation)), 0> type;
};

/*!
 * \brief Set the property for the material law by retrieving it from
 *        the spatial parameters.
 */
SET_PROP(BoxTwoP, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;

public:
    typedef typename SpatialParameters::MaterialLaw type;
};

/*!
 * \brief Set the property for the material parameters by extracting
 *        it from the material law.
 */
SET_PROP(BoxTwoP, MaterialLawParams)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLaw)) MaterialLaw;

public:
    typedef typename MaterialLaw::Params type;
};

SET_TYPE_PROP(BoxTwoP, FluidSystem, FluidSystem2P<TypeTag>);

SET_TYPE_PROP(BoxTwoP, FluidState, TwoPFluidState<TypeTag>);

// enable jacobian matrix recycling by default
SET_BOOL_PROP(BoxTwoP, EnableJacobianRecycling, true);
// enable partial reassembling by default
SET_BOOL_PROP(BoxTwoP, EnablePartialReassemble, true);

// \}
}

}

#endif
