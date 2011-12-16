/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
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
 * \ingroup Properties
 * \ingroup BoxProperties
 * \ingroup RichardsModel
 * \file
 *
 * \brief Contains the default definitions for the properties required
 *        by the Richards box model.
 */
#ifndef DUMUX_RICHARDS_PROPERTY_DEFAULTS_HH
#define DUMUX_RICHARDS_PROPERTY_DEFAULTS_HH

#include "richardsmodel.hh"
#include "richardsproblem.hh"
#include "richardsindices.hh"
#include "richardsfluxvariables.hh"
#include "richardsvolumevariables.hh"
#include "richardsproperties.hh"
#include "richardsnewtoncontroller.hh"

#include <dumux/material/components/nullcomponent.hh>
#include <dumux/material/fluidsystems/2pimmisciblefluidsystem.hh>

namespace Dumux
{
// \{

namespace Properties {
//////////////////////////////////////////////////////////////////
// Properties values
//////////////////////////////////////////////////////////////////
//! Number of equations required by the model
SET_INT_PROP(BoxRichards, NumEq, 1);
//! Number of fluid phases considered
SET_INT_PROP(BoxRichards, NumPhases, 2);

//! The local residual operator
SET_TYPE_PROP(BoxRichards,
              LocalResidual,
              RichardsLocalResidual<TypeTag>);

//! The global model used
SET_TYPE_PROP(BoxRichards, Model, RichardsModel<TypeTag>);

//! The class for the volume averaged quantities
SET_TYPE_PROP(BoxRichards, VolumeVariables, RichardsVolumeVariables<TypeTag>);

//! The class for the quantities required for the flux calculation
SET_TYPE_PROP(BoxRichards, FluxVariables, RichardsFluxVariables<TypeTag>);

//! The class of the newton controller
SET_TYPE_PROP(BoxRichards, NewtonController, RichardsNewtonController<TypeTag>);

//! The upwind weight for the mass conservation equations
SET_SCALAR_PROP(BoxRichards,
                MassUpwindWeight,
                1.0);

//! The class with all index definitions for the model
SET_TYPE_PROP(BoxRichards, RichardsIndices, Dumux::RichardsIndices);

/*!
 * \brief Set type of the parameter objects for the material law
 *
 * By default this is just retrieved from the material law.
 */
SET_PROP(BoxRichards, MaterialLawParams)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLaw)) MaterialLaw;

public:
    typedef typename MaterialLaw::Params type;
};

/*!
 * \brief The wetting phase used.
 *
 * By default we use the null-phase, i.e. this has to be defined by
 * the problem for the program to work. Please be aware that you
 * should be careful to use the Richards model in conjunction with
 * liquid non-wetting phases. This is only meaningful if the viscosity
 * of the liquid phase is _much_ lower than the viscosity of the
 * wetting phase.
 */
SET_PROP(BoxRichards, WettingPhase)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::NullComponent<Scalar> > type;
};

/*!
 * \brief The wetting phase used.
 *
 * By default we use the null-phase, i.e. this has to be defined by
 * the problem for the program to work. This doed not need to be
 * specified by the problem for the Richards model to work because the
 * Richards model does not conserve the non-wetting phase.
 */
SET_PROP(BoxRichards, NonWettingPhase)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::GasPhase<Scalar, Dumux::NullComponent<Scalar> > type;
};

/*!
 *\brief The fluid system used by the model.
 *
 * By default this uses the immiscible twophase fluid system. The
 * actual fluids used are specified using in the problem definition by
 * the WettingPhase and NonWettingPhase properties. Be aware that
 * using different fluid systems in conjunction with the Richards
 * model only makes very limited sense.
 */
SET_PROP(BoxRichards, FluidSystem)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(WettingPhase)) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NonWettingPhase)) NonWettingPhase;

public:
    typedef Dumux::FluidSystems::TwoPImmiscible<Scalar,
                                                WettingPhase,
                                                NonWettingPhase> type;
};

// \}
};

} // end namepace

#endif
