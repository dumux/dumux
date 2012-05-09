// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2010 by Andreas Lauser                               *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \ingroup TwoPBoxModel
 * \ingroup BoxProperties
 * \ingroup Properties
 * \file
 *
 * \brief Defines default values for the properties required by the
 *        twophase box model.
 */
#ifndef DUMUX_2PDFM_PROPERTY_DEFAULTS_HH
#define DUMUX_2PDFM_PROPERTY_DEFAULTS_HH

#include "2pdfmmodel.hh"
#include "2pdfmproblem.hh"
#include "2pdfmindices.hh"
#include "2pdfmfluxvariables.hh"
#include "2pdfmvolumevariables.hh"
#include "2pdfmproperties.hh"

#include <dumux/material/fluidsystems/gasphase.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/nullcomponent.hh>

#include <dumux/material/fluidsystems/2pimmisciblefluidsystem.hh>
#include <dumux/material/fluidstates/immisciblefluidstate.hh>

namespace Dumux
{
namespace Properties
{
//////////////////////////////////////////////////////////////////
// Property defaults
//////////////////////////////////////////////////////////////////
SET_INT_PROP(BoxTwoPDFM, NumEq, 2); //!< set the number of equations to 2
SET_INT_PROP(BoxTwoPDFM, NumPhases, 2); //!< The number of phases in the 2p model is 2

//! Set the default formulation to pWsN
SET_INT_PROP(BoxTwoPDFM,
             Formulation,
             TwoPCommonIndices::pwSn);

//! Use the 2p local jacobian operator for the 2p model
SET_TYPE_PROP(BoxTwoPDFM,
              LocalResidual,
              TwoPLocalResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(BoxTwoPDFM, Model, TwoPDFMModel<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(BoxTwoPDFM, VolumeVariables, TwoPDFMVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(BoxTwoPDFM, FluxVariables, TwoPDFMFluxVariables<TypeTag>);

//! the upwind weight for the mass conservation equations.
SET_SCALAR_PROP(BoxTwoPDFM, MassUpwindWeight, 1.0);

//! The indices required by the isothermal 2p model
SET_TYPE_PROP(BoxTwoPDFM,
              Indices,
              TwoPDFMIndices<GET_PROP_VALUE(TypeTag, Formulation), 0>);

//! DEPRECATED TwoPDFMIndices property
SET_TYPE_PROP(BoxTwoPDFM, TwoPDFMIndices, typename GET_PROP_TYPE(TypeTag, Indices));

//! DEPRECATED SpatialParameters property
SET_TYPE_PROP(BoxTwoPDFM, SpatialParameters, typename GET_PROP_TYPE(TypeTag, SpatialParams));

/*!
 * \brief Set the property for the material parameters by extracting
 *        it from the material law.
 */
SET_TYPE_PROP(BoxTwoPDFM,
              MaterialLawParams,
              typename GET_PROP_TYPE(TypeTag, MaterialLaw)::Params);

SET_PROP(BoxTwoPDFM, WettingPhase)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::NullComponent<Scalar> > type;
};

SET_PROP(BoxTwoPDFM, NonwettingPhase)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::NullComponent<Scalar> > type;
};

SET_PROP(BoxTwoPDFM, FluidSystem)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, NonwettingPhase) NonwettingPhase;

public:
    typedef Dumux::FluidSystems::TwoPImmiscible<Scalar,
                                                WettingPhase,
                                                NonwettingPhase> type;
};

SET_PROP(BoxTwoPDFM, FluidState)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
public:
    typedef ImmiscibleFluidState<Scalar, FluidSystem> type;
};

// disable velocity output by default
SET_BOOL_PROP(BoxTwoPDFM, EnableVelocityOutput, false);

}
//

}

#endif
