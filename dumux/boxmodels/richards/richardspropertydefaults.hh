// $Id$
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

//! The weight of the upwind vertex when looking at the mobility
SET_SCALAR_PROP(BoxRichards,
                MobilityUpwindAlpha,
                1.0);

//! The class with all index definitions for the model
SET_TYPE_PROP(BoxRichards, RichardsIndices, Dumux::RichardsIndices);

/*!
 * \brief The material law for capillary pressure and relative permeability
 *
 * By default this type is determined by retrieving it from the
 * spatial parameters.
 */
SET_PROP(BoxRichards, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;

public:
    typedef typename SpatialParameters::MaterialLaw type;
};

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
 *\brief The fluid system used by the model.
 *
 * By default this uses the immiscible twophase fluid system. The
 * actual fluids used are specified using in the problem definition by
 * the WettingPhase and NonwettingPhase properties. Be aware that
 * using different fluid systems in conjunction with the Richards
 * model only makes very limited sense.
 */
SET_TYPE_PROP(BoxRichards, FluidSystem, FluidSystem2P<TypeTag>);

/*!
 * \brief The non-wetting phase used.
 *
 * By default we use gaseous nitrogen as non-wetting phase. Please be
 * aware that you should be careful to use the Richards model in
 * conjunction with liquid non-wetting phases. This is only meaningful
 * if the viscosity of the liquid phase is _much_ lower than the
 * viscosity of the wetting phase.
 */
SET_PROP(BoxRichards, NonwettingPhase)
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef GasPhase<Scalar, N2<Scalar> > type;
};

//! The fluid state class
SET_TYPE_PROP(BoxRichards, FluidState, RichardsFluidState<TypeTag>);

// \}
};

} // end namepace

#endif
