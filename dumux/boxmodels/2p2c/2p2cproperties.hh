// $Id: 2p2cproperties.hh 3784 2010-06-24 13:43:57Z bernd $
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
 * \brief Defines the properties required for the 2p2c BOX model.
 */

#ifndef DUMUX_2P2C_PROPERTIES_HH
#define DUMUX_2P2C_PROPERTIES_HH

#include <dumux/boxmodels/common/boxproperties.hh>

namespace Dumux
{
/*!
 * \addtogroup TwoPTwoCModel
 */
// \{
namespace Properties
{
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the isothermal single phase problems
NEW_TYPE_TAG(BoxTwoPTwoC, INHERITS_FROM(BoxModel));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(NumPhases);   //!< Number of fluid phases in the system
NEW_PROP_TAG(NumComponents); //!< Number of fluid components in the system
NEW_PROP_TAG(TwoPTwoCIndices); //!< Enumerations for the 2p2c models
NEW_PROP_TAG(Formulation);   //!< The formulation of the model
NEW_PROP_TAG(SpatialParameters); //!< The type of the spatial parameters
NEW_PROP_TAG(FluidSystem); //!< Type of the multi-component relations

NEW_PROP_TAG(MaterialLaw);   //!< The material law which ought to be used (extracted from the spatial parameters)
NEW_PROP_TAG(MaterialLawParams); //!< The context material law (extracted from the spatial parameters)

NEW_PROP_TAG(EnableGravity); //!< Returns whether gravity is considered in the problem
NEW_PROP_TAG(MobilityUpwindAlpha); //!< The value of the upwind parameter for the mobility
}
}

namespace Dumux {

////////////////////////////////
// forward declarations
////////////////////////////////
template<class TypeTag>
class TwoPTwoCModel;

template<class TypeTag>
class TwoPTwoCLocalResidual;

template <class TypeTag>
class TwoPTwoCVolumeVariables;

template <class TypeTag>
class TwoPTwoCFluxVariables;

template <class TypeTag>
class TwoPTwoCNewtonController;

/*!
 * \brief Enumerates the formulations which the 2p2c model accepts.
 */
struct TwoPTwoCFormulation
{
    enum {
        plSg,
        pgSl
    };
};

/*!
 * \brief The indices for the isothermal TwoPTwoC model.
 *
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <class TypeTag,
          int formulation = TwoPTwoCFormulation::plSg,
          int PVOffset = 0>
class TwoPTwoCIndices
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
    // Phase indices
    static const int lPhaseIdx = FluidSystem::lPhaseIdx; //!< Index of the liquid phase
    static const int gPhaseIdx = FluidSystem::gPhaseIdx; //!< Index of the gas phase

    // Component indices
    static const int lCompIdx = 0; //!< Index of the liquid's primary component
    static const int gCompIdx = 1; //!< Index of the gas' primary component

    // present phases (-> 'pseudo' primary variable)
    static const int lPhaseOnly = 1; //!< Only the non-wetting phase is present
    static const int gPhaseOnly = 0; //!< Only the wetting phase is present
    static const int bothPhases = 2; //!< Both phases are present

    // Primary variable indices
    static const int pressureIdx = PVOffset + 0; //!< Index for wetting/non-wetting phase pressure (depending on formulation) in a solution vector
    static const int switchIdx = PVOffset + 1; //!< Index of the either the saturation or the mass fraction of the non-wetting/wetting phase

    static const int plIdx = pressureIdx; //!< Index for liquid phase pressure in a solution vector
    static const int SgOrXIdx = switchIdx; //!< Index of the either the saturation of the gas phase or the mass fraction secondary component in the only phase

    // equation indices
    static const int contiLEqIdx = PVOffset + 0; //!< Index of the mass conservation equation for the liquid's primary component
    static const int contiGEqIdx = PVOffset + 1; //!< Index of the mass conservation equation for the gas' primary component
};

/*!
 * \brief The indices for the isothermal TwoPTwoC model in the pg-Sl
 *        formulation.
 *
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <class TypeTag, int PVOffset>
class TwoPTwoCIndices<TypeTag, TwoPTwoCFormulation::pgSl, PVOffset>
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
    // Phase indices
    static const int lPhaseIdx = FluidSystem::lPhaseIdx; //!< Index of the liquid phase
    static const int gPhaseIdx = FluidSystem::gPhaseIdx; //!< Index of the gas phase

    // Component indices
    static const int lCompIdx = 0; //!< Index of the liquid's primary component
    static const int gCompIdx = 1; //!< Index of the gas' primary component

    // present phases (-> 'pseudo' primary variable)
    static const int lPhaseOnly = 1; //!< Only the non-wetting phase is present
    static const int gPhaseOnly = 2; //!< Only the wetting phase is present
    static const int bothPhases = 3; //!< Both phases are present

    // Primary variable indices
    static const int pressureIdx = PVOffset + 0; //!< Index for wetting/non-wetting phase pressure (depending on formulation) in a solution vector
    static const int switchIdx = PVOffset + 1; //!< Index of the either the saturation or the mass fraction of the non-wetting/wetting phase

    static const int pgIdx = pressureIdx; //!< Index for gas phase pressure in a solution vector
    static const int SlOrXIdx = switchIdx; //!< Index of the either the saturation of the liquid phase or the mass fraction secondary component in the only phase

    // Equation indices
    static const int contiLEqIdx = PVOffset + 1; //!< Index of the mass conservation equation for the liquid's primary component
    static const int contiGEqIdx = PVOffset + 0; //!< Index of the mass conservation equation for the gas' primary component
};



namespace Properties {

//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

/*!
 * \brief Set the property for the number of components.
 *
 * We just forward the number from the fluid system and use an static
 * assert to make sure it is 2.
 */
SET_PROP(BoxTwoPTwoC, NumComponents)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
    static const int value = FluidSystem::numComponents;

    static_assert(value == 2,
                  "Only fluid systems with 2 components are supported by the 2p-2c model!");
};

/*!
 * \brief Set the property for the number of fluid phases.
 *
 * We just forward the number from the fluid system and use an static
 * assert to make sure it is 2.
 */
SET_PROP(BoxTwoPTwoC, NumPhases)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
    static const int value = FluidSystem::numPhases;
    static_assert(value == 2,
                  "Only fluid systems with 2 phases are supported by the 2p-2c model!");
};

SET_INT_PROP(BoxTwoPTwoC, NumEq, 2); //!< set the number of equations to 2

//! Set the default formulation to pl-Sg
SET_INT_PROP(BoxTwoPTwoC,
             Formulation,
             TwoPTwoCFormulation::plSg);

/*!
 * \brief Set the property for the material law by retrieving it from
 *        the spatial parameters.
 */
SET_PROP(BoxTwoPTwoC, MaterialLaw)
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
SET_PROP(BoxTwoPTwoC, MaterialLawParams)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLaw)) MaterialLaw;

public:
    typedef typename MaterialLaw::Params type;
};

//! Use the 2p2c local jacobian operator for the 2p2c model
SET_TYPE_PROP(BoxTwoPTwoC,
              LocalResidual,
              TwoPTwoCLocalResidual<TypeTag>);

//! Use the 2p2c specific newton controller for the 2p2c model
SET_TYPE_PROP(BoxTwoPTwoC, NewtonController, TwoPTwoCNewtonController<TypeTag>);

//! the Model property
SET_TYPE_PROP(BoxTwoPTwoC, Model, TwoPTwoCModel<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(BoxTwoPTwoC, VolumeVariables, TwoPTwoCVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(BoxTwoPTwoC, FluxVariables, TwoPTwoCFluxVariables<TypeTag>);

//! the upwind factor for the mobility.
SET_SCALAR_PROP(BoxTwoPTwoC, MobilityUpwindAlpha, 1.0);

//! The indices required by the isothermal 2p2c model
SET_PROP(BoxTwoPTwoC,
         TwoPTwoCIndices)
{ private:
    enum { Formulation = GET_PROP_VALUE(TypeTag, PTAG(Formulation)) };
public:
    typedef TwoPTwoCIndices<TypeTag, Formulation, 0> type;
};

// enable jacobian matrix recycling by default
SET_BOOL_PROP(BoxTwoPTwoC, EnableJacobianRecycling, true);
// enable partial reassembling by default
SET_BOOL_PROP(BoxTwoPTwoC, EnablePartialReassemble, true);
// enable time-step ramp up by default
SET_BOOL_PROP(BoxTwoPTwoC, EnableTimeStepRampUp, true);

// \}
}

}

#endif
