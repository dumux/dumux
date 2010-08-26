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
 * \brief Contains the default definitions for the properties for the
 *        Richards BOX model.
 */
#ifndef DUMUX_RICHARDS_PROPERTY_DEFAULTS_HH
#define DUMUX_RICHARDS_PROPERTY_DEFAULTS_HH

#include "richardsmodel.hh"
#include "richardsproblem.hh"
#include "richardsindices.hh"
#include "richardsfluxvariables.hh"
#include "richardsvolumevariables.hh"
#include "richardsfluidstate.hh"
#include "richardsproperties.hh"
#include "richardsnewtoncontroller.hh"

#include <dumux/material/fluidsystems/2p_system.hh>
#include <dumux/material/components/n2.hh>

namespace Dumux
{
/*!
 * \addtogroup RichardsModel
 */
// \{

namespace Properties {
//////////////////////////////////////////////////////////////////
// Properties values
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

//! the NewtonController property
SET_TYPE_PROP(BoxRichards, NewtonController, RichardsNewtonController<TypeTag>);

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

SET_TYPE_PROP(BoxRichards, FluidSystem, FluidSystem2P<TypeTag>);
SET_PROP(BoxRichards, NonwettingPhase)
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef GasPhase<Scalar, N2<Scalar>> type;
};

SET_TYPE_PROP(BoxRichards, FluidState, RichardsFluidState<TypeTag>);

// \}
};

} // end namepace

#endif
