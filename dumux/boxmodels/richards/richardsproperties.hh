// $Id: richardsproperties.hh 3840 2010-07-15 10:14:15Z bernd $
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
#ifndef DUMUX_RICHARDS_PROPERTIES_HH
#define DUMUX_RICHARDS_PROPERTIES_HH

namespace Dumux
{
/*!
 * \addtogroup RichardsModel
 */
// \{
////////////////////////////////
// forward declarations
////////////////////////////////
template<class TypeTag>
class RichardsModel;

template<class TypeTag>
class RichardsLocalResidual;

template <class TypeTag>
class RichardsVolumeVariables;

template <class TypeTag>
class RichardsFluxVariables;

/*!
 * \brief Indices for the single phase model.
 */
struct RichardsIndices
{
    //! index of the wetting phase pressure
    static const int pW = 0;
};

///////////////////////////////////////////////////////////////////////////
// properties for the isothermal richards model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for problems discretized using the isothermal
//! richards model
NEW_TYPE_TAG(BoxRichards, INHERITS_FROM(BoxModel));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(NumPhases);   //!< Number of fluid phases in the system
NEW_PROP_TAG(RichardsIndices); //!< Enumerations for the richards models
NEW_PROP_TAG(SpatialParameters); //!< The type of the spatial parameters object
NEW_PROP_TAG(MaterialLaw);   //!< The material law which ought to be used (extracted from the spatial parameters)
NEW_PROP_TAG(MaterialLawParams); //!< The context material law (extracted from the spatial parameters)
NEW_PROP_TAG(WettingPhase); //!< The wetting phase for the richards model
NEW_PROP_TAG(EnableGravity); //!< Returns whether gravity is considered in the problem
NEW_PROP_TAG(MobilityUpwindAlpha); //!< The value of the upwind parameter for the mobility

//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////
SET_INT_PROP(BoxRichards, NumEq, 1);
SET_INT_PROP(BoxRichards, NumPhases, 2);

//! Use the 2p local jacobian operator for the 2p model
SET_TYPE_PROP(BoxRichards,
              LocalResidual,
              RichardsLocalResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(BoxRichards, Model, RichardsModel<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(BoxRichards, VolumeVariables, RichardsVolumeVariables<TypeTag>);




//! the FluxVariables property
SET_TYPE_PROP(BoxRichards, FluxVariables, RichardsFluxVariables<TypeTag>);

//! the weight of the upwind vertex for the mobility
SET_SCALAR_PROP(BoxRichards,
                MobilityUpwindAlpha,
                1.0);

//! The indices required by the isothermal single-phase model
SET_TYPE_PROP(BoxRichards, RichardsIndices, Dumux::RichardsIndices);

/*!
 * \brief Set the property for the material law by retrieving it from
 *        the spatial parameters.
 */
SET_PROP(BoxRichards, MaterialLaw)
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
SET_PROP(BoxRichards, MaterialLawParams)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLaw)) MaterialLaw;

public:
    typedef typename MaterialLaw::Params type;
};

// \}
};

} // end namepace

#endif
