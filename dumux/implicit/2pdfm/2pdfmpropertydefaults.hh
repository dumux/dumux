/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Defines default values for the properties required by the
 *        2pDFM fully implicit model.
 *
 * \ingroup Properties
 * \ingroup TwoPDFMModel
 * \ingroup ImplicitProperties
 */
#ifndef DUMUX_MODELS_2PDFM_PROPERTY_DEFAULTS_HH
#define DUMUX_MODELS_2PDFM_PROPERTY_DEFAULTS_HH

#include <dumux/porousmediumflow/implicit/darcyfluxvariables.hh>
#include <dumux/material/components/nullcomponent.hh>
#include <dumux/material/fluidstates/immisciblefluidstate.hh>
#include <dumux/material/fluidsystems/2pimmisciblefluidsystem.hh>
#include <dumux/material/fluidsystems/gasphase.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/spatialparams/implicitspatialparams.hh>

#include "2pdfmmodel.hh"
#include "2pdfmproblem.hh"
#include "2pdfmindices.hh"
#include "2pdfmfluxvariables.hh"
#include "2pdfmvolumevariables.hh"
#include "2pdfmproperties.hh"
#include "2pdfmlocalresidual.hh"

namespace Dumux
{
namespace Properties
{
//////////////////////////////////////////////////////////////////
// Property defaults
//////////////////////////////////////////////////////////////////
SET_INT_PROP(TwoPDFM, NumEq, 2); //!< set the number of equations to 2
SET_INT_PROP(TwoPDFM, NumPhases, 2); //!< The number of phases in the 2p model is 2

//! Set the default formulation to pWsN
SET_INT_PROP(TwoPDFM,
             Formulation,
             TwoPFormulation::pwsn);

//! Use the 2p local jacobian operator for the 2p model
SET_TYPE_PROP(TwoPDFM,
              LocalResidual,
              TwoPDFMLocalResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(TwoPDFM, Model, TwoPDFMModel<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(TwoPDFM, VolumeVariables, TwoPDFMVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(TwoPDFM, FluxVariables, TwoPDFMFluxVariables<TypeTag>);

//! the upwind weight for the mass conservation equations.
SET_SCALAR_PROP(TwoPDFM, ImplicitMassUpwindWeight, 1.0);

//! weight for the upwind mobility in the velocity calculation
SET_SCALAR_PROP(TwoPDFM, ImplicitMobilityUpwindWeight, 1.0);

//! The indices required by the isothermal 2pDFM model
SET_PROP(TwoPDFM, Indices)
{ private:
    enum { Formulation = GET_PROP_VALUE(TypeTag, Formulation) };
 public:
    typedef TwoPDFMIndices<TypeTag, Formulation, 0> type;
};


//! The spatial parameters to be employed.
//! Use ImplicitSpatialParams by default.
SET_TYPE_PROP(TwoPDFM, SpatialParams, ImplicitSpatialParams<TypeTag>);

/*!
 * \brief Set the property for the material parameters by extracting
 *        it from the material law.
 */
SET_TYPE_PROP(TwoPDFM,
              MaterialLawParams,
              typename GET_PROP_TYPE(TypeTag, MaterialLaw)::Params);

SET_PROP(TwoPDFM, WettingPhase)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::NullComponent<Scalar> > type;
};

SET_PROP(TwoPDFM, NonwettingPhase)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::NullComponent<Scalar> > type;
};

SET_PROP(TwoPDFM, FluidSystem)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, NonwettingPhase) NonwettingPhase;

public:
    typedef Dumux::FluidSystems::TwoPImmiscible<Scalar,
                                                WettingPhase,
                                                NonwettingPhase> type;
};

SET_PROP(TwoPDFM, FluidState)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
public:
    typedef ImmiscibleFluidState<Scalar, FluidSystem> type;
};

// disable velocity output by default
SET_BOOL_PROP(TwoPDFM, VtkAddVelocity, false);

// enable gravity by default
SET_BOOL_PROP(TwoPDFM, ProblemEnableGravity, true);
} // end namespace Properties
} // end namespace Dumux

#endif // DUMUX_MODELS_2PDFM_PROPERTY_DEFAULTS_HH
